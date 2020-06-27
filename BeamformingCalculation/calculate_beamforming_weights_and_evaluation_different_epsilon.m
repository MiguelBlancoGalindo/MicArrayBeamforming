function [BFEvalAllArrays, BFWeightsAllArrays, ConfigAllArrays] = calculate_beamforming_weights_and_evaluation_different_epsilon(ConfigAllArrays,ArrayManAllArrays,SourceDirection,Interferers,epsilon)
% function [BFEvalAllArrays, BFWeightsAllArrays, ConfigAllArrays] = calculate_beamforming_weights_and_evaluation_different_epsilon(ConfigAllArrays,ArrayManAllArrays,SourceDirection,Interferers,epsilon)
% function to calculate and evaluate the beamforming performance of several
% methods. 
%
% input arguments: 
%   ConfigAllArrays: configuration struct containing all the settings.
%   ArrayManAllArrays: struct containing the array manifold transfer functions.
%   SourceDirection: L x 2 matrix with elevation and azimuth angles of L 
%       look directions.
%   InterfererDirection: L x I*2 matrix with pairs of elevation and azimuth 
%       angles of I interferer angles for each L look direction.
%   epsilon: regularisation parameter vector.
%
% output arguments:
%   BFEvalAllArrays: evaluation struct.
%   BFWeightsAllArrays: beamforming weights struct.
%   ConfigAllArrays: configuration struct containing all the settings.

ne = length(epsilon);
nArrayGeo = size(ConfigAllArrays,1);
nM = size(ConfigAllArrays,2);
nr = size(ConfigAllArrays,3);
BFEvalAllArrays = repmat(struct([]),nArrayGeo,nM,nr,ne);
BFWeightsAllArrays = repmat(struct([]),nArrayGeo,nM,nr,ne);
M = zeros(nM,1);
targain = 1;  
intgain = 1; 
diffgain = 0;

for iM=1:nM
    for ir=1:nr
        for iarray_geo=1:nArrayGeo
            for ie=1:ne
            ConfigAllArrays{iarray_geo,iM,ir,ie}=configFiltersEpsilon(ConfigAllArrays{iarray_geo,iM,ir},SourceDirection,Interferers,epsilon(ie));
            end
        end
    end
    M(iM) = size(ConfigAllArrays{iarray_geo,iM,ir,ie}.Array.MicPos,1);
end
%creating string of beamforming method to use when saving the results
nmethods = length(ConfigAllArrays{1,1,1}.Filter.Method);    %number of bf methods
%determining whether Rxx needs to be calculating depending on whether mvdr
%or lcmv are to be computed
getRxx=0;
Rxx=[];
SetupAudio=[];
for imethod=1:nmethods
    if strcmp(ConfigAllArrays{1,1,1,1}.Filter.Method{imethod},'mvdr') || ...
        strcmp(ConfigAllArrays{1,1,1,1}.Filter.Method{imethod},'lcmv')
        getRxx=1;
        Rxx = cell(size(ConfigAllArrays{1,1,1,1}.Scene.Targets.Angle,1),1);
    end
end

%calculating indices of targets and interferers
%needs to be wrapped around -180 and 180 degrees. Not done yet
tarind = zeros(size(ConfigAllArrays{1,1,1,1}.Scene.Targets.Angle,1),size(ConfigAllArrays{1,1,1,1}.Scene.Targets.Angle,2));
intind = zeros(size(ConfigAllArrays{1,1,1,1}.Scene.Interferers.Angle,1),size(ConfigAllArrays{1,1,1,1}.Scene.Interferers.Angle,2));
%same gain for all target and interferers
intgain = repmat(intgain,size(intind,2)/2,1);

for itar=1:size(ConfigAllArrays{1,1,1,1}.Scene.Targets.Angle,1)
    [~,tarind(itar,1)]=min(abs(rad2deg(ArrayManAllArrays{1,1,1,1}.Setup.LookAng.Az) - ConfigAllArrays{1,1,1,1}.Scene.Targets.Angle(itar,1)));
    if numel(tarind)>size(tarind,1)
        [~,tarind(itar,2)]=min(abs(rad2deg(ArrayManAllArrays{1,1,1,1}.Setup.LookAng.El) - ConfigAllArrays{1,1,1,1}.Scene.Targets.Angle(itar,2)));
    end
end
for itar=1:size(ConfigAllArrays{1,1,1,1}.Scene.Interferers.Angle,1)
    for iint=1:size(ConfigAllArrays{1,1,1,1}.Scene.Interferers.Angle,2)/2
    [~,intind(itar,iint*2-1)]=min(abs(rad2deg(ArrayManAllArrays{1,1,1}.Setup.LookAng.Az) - ConfigAllArrays{1,1,1,1}.Scene.Interferers.Angle(itar,iint*2-1)));
        if numel(intind)>size(intind,1)
            [~,intind(itar,iint*2)]=min(abs(rad2deg(ArrayManAllArrays{1,1,1}.Setup.LookAng.El) - ConfigAllArrays{1,1,1,1}.Scene.Interferers.Angle(itar,iint*2)));
        end
    end
end

%calculating Rxx if needed, the weights of the beamforming method and the
%metrics evaluating its performance
nInts = size(ConfigAllArrays{1,1,1,1}.Scene.Interferers.Angle,2)/2;
nTars = size(ConfigAllArrays{1,1,1,1}.Scene.Targets.Angle,1);
freqs = size(ArrayManAllArrays{1,1,1}.Setup.a,4);
nAngles = length(ArrayManAllArrays{1,1,1}.Setup.LookAng.Az)*length(ArrayManAllArrays{1,1,1}.Setup.LookAng.El);
if strcmpi(ConfigAllArrays{1,1,1,1}.Filter.RegMethod,'external')
    RegParamStr = 'Epsilon';
elseif strcmpi(ConfigAllArrays{1,1,1,1}.Filter.RegMethod,'wnglimit')
    RegParamStr = 'WNGmin';
end
for iM=1:nM
    for ir=1:nr
        for iarray_geo=1:nArrayGeo
            % Calculating cross-correlation matrix of input signal Rxx from target and
            % interferer plane waves (rather than audio signals) for MVDR
            if getRxx
                Rii = zeros(M(iM),M(iM),freqs);
                Rxx = cell(size(ConfigAllArrays{1,1,1,1}.Scene.Targets.Angle,1),1);
                Rnn = cell(size(ConfigAllArrays{1,1,1,1}.Scene.Targets.Angle,1),1);
                xtar = zeros(M(iM),nTars);
                Rintint = zeros(nInts,1);
                for iTar=1:nTars
                    Rxx{iTar} = zeros(M(iM),M(iM),freqs);
                    A_s = zeros(M(iM),nInts+1,freqs);
                    for f=1:freqs
                        xint = zeros(M(iM),nInts);
                        ilook=1;
                        for iInt=1:nInts
                            intarrayman=squeeze(ArrayManAllArrays{iarray_geo,iM,ir}.Setup.a(:,intind(iTar,iInt*2),intind(iTar,iInt*2-1),f));
                            xint(:,iInt) = intarrayman*intgain(iInt);
                            Rintint(iInt) = intgain*intgain';
                            ilook = ilook+1;
                            A_s(:,ilook,f) = intarrayman;
                        end
                        Rtartar = targain*targain';
                        xint = sum(xint,2);
                        tararrayman=squeeze(ArrayManAllArrays{iarray_geo,iM,ir}.Setup.a(:,tarind(iTar,2),tarind(iTar,1),f));
                        xtar(:,iTar) = tararrayman*targain;
                        x = xtar(:,iTar) + xint;
                        Rii(:,:,f)=xint*xint';
                        Rxx{iTar}(:,:,f)=x*x';
                        A_s(:,1,f) = tararrayman;
                        
                        %Rxx{iTar}(:,:,f)=Rii(:,:,f);
                    end
                    Rss = diag([Rtartar; Rintint]);
                    Rvv = diffgain.*squeeze(ArrayManAllArrays{iarray_geo,iM,ir}.Setup.Rvv);
                    Rnn{iTar} = Rvv + Rii;
                    %Rxx{iTar} = Rxx{iTar} + Rvv;
                    for f=1:freqs
                        A_sk = A_s(:,:,f);
                        Rxx{iTar}(:,:,f) = A_sk*Rss*A_sk';
                    end
                    Rxx{iTar} = Rxx{iTar} + Rvv;
                end
            end    
            
            for ie=1:ne
            [BFWeights,BFEval]=RunBFTests(ConfigAllArrays{iarray_geo,iM,ir,ie},ArrayManAllArrays{iarray_geo,iM,ir},SetupAudio,Rxx); % use [] for 3rd argument if no audio provided (can't use MDVR or LCMV in this case)
            if ~isfield(ConfigAllArrays{1,1,1,1}.ArrayMan, 'Ntests') || ConfigAllArrays{1,1,1,1}.ArrayMan.Ntests==1
                BFEvalAllArrays{iarray_geo,iM,ir,ie} = BFEval;
                BFWeightsAllArrays{iarray_geo,iM,ir,ie} = BFWeights;
            end
            disp(['Weights and evaluation calculated for regularization ' num2str(ie) ' of ' num2str(ne)]);
            disp([num2str(((iarray_geo-1)*ne+ie)/(nArrayGeo*ne)*100) '%']);
            Config = ConfigAllArrays{iarray_geo,iM,ir,ie};
            if Config.Temp.DoSave
                save([Config.Temp.ResultsFolder 'Beamforming_N5_' Config.Array.Geo '_' RegParamStr num2str(ie) '_log20Hz-20kHz_' Config.Temp.SimName '_' num2str(Config.ArrayMan.Ntests) 'tests.mat'], 'BFEval','BFWeights','Config','-v7.3');
            end
            end
        end

    end
end

end