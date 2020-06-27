function [BFEvalAllArrays, BFWeightsAllArrays, ConfigAllArrays] = calculate_beamforming_weights_and_evaluation_different_arrays(ConfigAllArrays,ArrayManAllArrays)
% function [BFEvalAllArrays, BFWeightsAllArrays, ConfigAllArrays] = calculate_beamforming_weights_and_evaluation_different_arrays(ConfigAllArrays,ArrayManAllArrays)
% function to calculate and evaluate the beamforming performance of several
% methods. 
%
% input arguments: 
%   ConfigAllArrays: configuration struct containing all the settings.
%   ArrayManAllArrays: struct containing the array manifold transfer functions.
% output arguments:
%   BFEvalAllArrays: evaluation struct.
%   BFWeightsAllArrays: beamforming weights struct.
%   ConfigAllArrays: configuration struct containing all the settings.

nArrayGeo = size(ConfigAllArrays,1);
nM = size(ConfigAllArrays,2);
nr = size(ConfigAllArrays,3);
BFEvalAllArrays = repmat(struct([]),nArrayGeo,nM,nr);
BFWeightsAllArrays = repmat(struct([]),nArrayGeo,nM,nr);
M = zeros(nM,1);

%creating string of beamforming method to use when saving the results
nmethods = length(ConfigAllArrays{1,1,1}.Filter.Method);    %number of bf methods
%determining whether Rxx needs to be calculating depending on whether mvdr
%or lcmv are to be computed
getRxx=0;
Rxx=[];
SetupAudio=[];
tarAngle = ConfigAllArrays{1,1,1}.Scene.Targets.Angle;
intAngle = ConfigAllArrays{1,1,1}.Scene.Interferers.Angle;
nTars = size(ConfigAllArrays{1,1,1}.Scene.Targets.Angle,1);
nInts = size(ConfigAllArrays{1,1,1}.Scene.Interferers.Angle,2)/2;

for imethod=1:nmethods
    if strcmp(ConfigAllArrays{1,1,1}.Filter.Method{imethod},'mvdr') || ...
        strcmp(ConfigAllArrays{1,1,1}.Filter.Method{imethod},'lcmv')
        getRxx=1;
        Rxx = cell(nTars,1);
    end
end

tarInd = zeros(nTars,1);
intInd = zeros(nTars,nInts);

targetAz = tarAngle(:,2);
targetEl = tarAngle(:,1);

%find indices of target and interferers
if ~isfield(ConfigAllArrays{1,1,1}.Scene.Targets,'AngleInd')
    for iTar=1:nTars
        [~,~,tarInd(iTar)] = getNearestAngle(ConfigAllArrays{1},targetEl(iTar),targetAz(iTar));
    end
end

if ~isfield(ConfigAllArrays{1,1,1}.Scene.Interferers,'AngleInd')
for iTar=1:nTars
    for iInt=1:nInts
        [~,~,intInd(iTar,iInt)] = getNearestAngle(ConfigAllArrays{1},intAngle(iInt*2-1),intAngle(iInt*2));
    end
end
end

%calculating Rxx if needed, the weights of the beamforming method and the
%metrics evaluating its performance
nfreqs = size(ArrayManAllArrays{1,1,1}.Setup.a,3);
for iM=1:nM
    M(iM) = ConfigAllArrays{1,iM,1}.Array.M;
    for ir=1:nr
        for iarray=1:nArrayGeo
            % Calculating cross-correlation matrix of input signal Rxx from target and
            % interferer plane waves (rather than audio signals) for MVDR
            if getRxx
                Rii = zeros(M(iM),M(iM),nfreqs);
                Rxx = cell(size(ConfigAllArrays{1,1,1}.Scene.Targets.Angle,1),1);
                Rnn = cell(size(ConfigAllArrays{1,1,1}.Scene.Targets.Angle,1),1);
                xtar = zeros(M(iM),nTars);
                Rintint = zeros(nInts,1);
                for iTar=1:nTars
                    Rxx{iTar} = zeros(M(iM),M(iM),nfreqs);
                    A_s = zeros(M(iM),nInts+1,nfreqs);
                    for f=1:nfreqs
                        xint = zeros(M(iM),nInts);
                        ilook=1;
                        for iInt=1:nInts
                            intarrayman=squeeze(ArrayManAllArrays{iarray,iM,ir}.Setup.a(:,intInd(iTar,iInt),f));
                            xint(:,iInt) = intarrayman*ConfigAllArrays{iarray,iM,ir}.Scene.Interferers.Gain(iInt);
                            Rintint(iInt) = ConfigAllArrays{iarray,iM,ir}.Scene.Interferers.Gain*ConfigAllArrays{iarray,iM,ir}.Scene.Interferers.Gain';
                            ilook = ilook+1;
                            A_s(:,ilook,f) = intarrayman;
                        end
                        Rtartar = ConfigAllArrays{iarray,iM,ir}.Scene.Targets.Gain*ConfigAllArrays{iarray,iM,ir}.Scene.Targets.Gain';
                        xint = sum(xint,2);
                        tararrayman=squeeze(ArrayManAllArrays{iarray,iM,ir}.Setup.a(:,tarInd(iTar),f));
                        xtar(:,iTar) = tararrayman*ConfigAllArrays{iarray,iM,ir}.Scene.Targets.Gain(iTar);
                        x = xtar(:,iTar) + xint;
                        Rii(:,:,f)=xint*xint';
                        Rxx{iTar}(:,:,f)=x*x';
                        A_s(:,1,f) = tararrayman;
                        
                        %Rxx{iTar}(:,:,f)=Rii(:,:,f);
                    end
                    Rss = diag([Rtartar; Rintint]);
                    diffGain = sqrt(ConfigAllArrays{iarray,iM,ir}.Scene.Diffuse.Coeff(iTar).*(ConfigAllArrays{iarray,iM,ir}.Scene.Targets.Gain(iTar).^2+sum(ConfigAllArrays{iarray,iM,ir}.Scene.Targets.Gain(iTar,:).^2)));
                    Rvv = diffGain.*squeeze(ArrayManAllArrays{iarray,iM,ir}.Setup.Rvv);
                    Rnn{iTar} = Rvv + Rii;
                    %Rxx{iTar} = Rxx{iTar} + Rvv;
                    for f=1:nfreqs
                        A_sk = A_s(:,:,f);
                        Rxx{iTar}(:,:,f) = A_sk*Rss*A_sk';
                    end
                    Rxx{iTar} = Rxx{iTar} + Rvv;
                end
            end    
            [BFWeights,BFEval]=RunBFTests(ConfigAllArrays{iarray,iM,ir},ArrayManAllArrays{iarray,iM,ir},SetupAudio,Rxx); % use [] for 3rd argument if no audio provided (can't use MDVR or LCMV in this case)
%             Config.Filter.CurrTarget=Config.Scene.Targets(1);
%             Config.Filter.CurrInts=Config.Scene.Interferers(:,1);
%             BFWeights=BeamformerWeights([],Config,ArrayMan,Config.Filter.Method{1});
%             BFEval=BFWeightEvaluation(BFWeights,Config,ArrayManAllArrays{1,m,n}.Monitor);
            BFEvalAllArrays{iarray,iM,ir} = BFEval;
            BFWeightsAllArrays{iarray,iM,ir} = BFWeights;
            
        end

    end
end


%save the beamforming methods
%save(['Results/Beamforming_el_fullbw_360SourceDir' methodStr '_fixed_' ArrayCriteria '_M' num2str(M) '_r' num2str(rmm) 'mm_' Config.ArrayMan.LookDistance '.mat'], 'BFEvalAllArrays','ConfigAllArrays','-v7.3');
end