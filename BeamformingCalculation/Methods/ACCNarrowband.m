function [ ACCWeights ] = ACCNarrowband(Config,ArrayMan)
% function [ ACCWeights ] = ACCNarrowband(Config,ArrayMan)
%
% Calculate narrowband filter weights based on acoustic contrast
% optimization on the parameters specified in the Config struct.
%
% input arguments: 
%   Config: configuration struct containing all the settings.
%   ArrayMan: struct containing the array manifold transfer functions.
%
% output arguments: 
%   ACCWeights: beamforming weights struct. 

doDeltaPhi = false;
if isfield(Config.Filter,'Acc') && isfield(Config.Filter.Acc,'PassBand') && ...
        isfield(Config.Filter.Acc,'StopBand')
    if isfield(Config.Filter.Acc.PassBand,'DeltaPhi') && ...
        isfield(Config.Filter.Acc.StopBand,'DeltaPhi')
        DeltaPhiPB = Config.Filter.Acc.PassBand.DeltaPhi;
        DeltaPhiSB = Config.Filter.Acc.StopBand.DeltaPhi;
        doDeltaPhi = true;
    elseif isfield(Config.Filter.Acc.PassBand,'Phi') && ...
        isfield(Config.Filter.Acc.StopBand,'Phi')
        PhiPB = Config.Filter.Acc.PassBand.Phi;
        PhiSB = Config.Filter.Acc.StopBand.Phi;
    else
        error('ACC setup parameters missing');
    end
else
    % added params for ACC-SDBF for PC book chapter figures
    warning('ACC hard-coded parameters correspond to PC book chapter settings');
    DeltaPhiPB=3;
    DeltaPhiSB=6;
end
if ~isfield(Config.Filter,'Acc') || ~isfield(Config.Filter.Acc,'Mode')
    Config.Filter.Acc.Mode='superdir'; % 'notch'
end

nFreqs=length(ArrayMan.f);
nMics=length(Config.Array.MicPos);
w=zeros(nMics,nFreqs);
I=eye(nMics);

targetAngle=Config.Filter.CurrTarget;
intAngles=Config.Filter.CurrInts;
nints=length(intAngles)/2;
[~,~,targetind] = getNearestAngle(Config,targetAngle(1),targetAngle(2));
if strcmp(Config.ArrayMan.SteerSpace,'2Daz')  
    if doDeltaPhi
        [~,~,acctargetind1] = getNearestAngle(Config,targetAngle(1),targetAngle(2)-DeltaPhiPB);
        [~,~,acctargetind2] = getNearestAngle(Config,targetAngle(1),targetAngle(2)+DeltaPhiPB);
        acctargetinds=acctargetind1:acctargetind2;
    else
        acctargetinds = zeros(length(PhiPB),1);
        for iPhi = 1:length(PhiPB)
            [~,~,acctargetinds(iPhi)] = getNearestAngle(Config,targetAngle(1),PhiPB(iPhi));
        end
    end

intinds = zeros(nints,1);
if strcmp(Config.Filter.Acc.Mode,'superdir')
    if doDeltaPhi
        [~,~,notintind1] = getNearestAngle(Config,targetAngle(1),targetAngle(2)-DeltaPhiSB);
        [~,~,notintind2] = getNearestAngle(Config,targetAngle(1),targetAngle(2)+DeltaPhiSB);
        notintinds=notintind1+1:notintind2-1;

        accintinds=1:length(Config.ArrayMan.LookAng.Az);
        accintinds(notintinds)=[];
        accintinds(accintinds<1)=[];
    else
        accintinds = zeros(length(PhiSB),1);
        for iPhi = 1:length(PhiSB)
            [~,~,accintinds(iPhi)] = getNearestAngle(Config,targetAngle(1),PhiSB(iPhi));
        end
    end
    
elseif strcmp(Config.Filter.Acc.Mode,'notch')
    
    for iints=1:nints
        [~,intinds(iints)]=min(abs(deg2rad(intAngles(iints))-ArrayMan.LookAng.Az));
        %[~,~,intinds(iints)] = getNearestAngle(Config,intAngles(kint*2-1),intAngles(kint*2));
        intindslong(iints,:)=intinds(iints)-2:intinds(iints)+2;
    end
    accintinds=unique(reshape(intindslong,1,numel(intindslong)));
    accintinds(accintinds<1)=[];
    
end

else
    error('ACC not implemented for sound fields other than cylindrical');
end

if length(Config.Filter.RegParam)==nFreqs
    if strcmp(Config.Filter.RegMethod,'wnglimit')
        WNGConst=Config.Filter.RegParam;
    elseif strcmp(Config.Filter.RegMethod,'external')
        epsilon=Config.Filter.RegParam;
    end
elseif length(Config.Filter.RegParam)==1
    if strcmp(Config.Filter.RegMethod,'wnglimit')
        WNGConst=repmat(Config.Filter.RegParam,nFreqs,1);
    elseif strcmp(Config.Filter.RegMethod,'external')
        epsilon=repmat(Config.Filter.RegParam,nFreqs,1);
    end
end

w(:,1)=zeros(nMics,1);
testepsilons=logspace(-6,3,200);   
wngtest=zeros(length(testepsilons),nFreqs);
for ifreq=1:nFreqs
    
    % Array manifold for that frequency and direction
    a_targets=ArrayMan.a(:,acctargetinds,ifreq);
    a_target=ArrayMan.a(:,targetind,ifreq);
    a_ints=ArrayMan.a(:,accintinds,ifreq);
    
    if all(all(isnan(a_target))); a_target = []; end
    if all(all(isnan(a_targets))); a_targets = []; end
    if all(all(isnan(a_ints))); a_ints = []; end
    
    % This is where we will investigate different definitions of Rb and Rd
    Rb=a_targets*a_targets';
    Rd=a_ints*a_ints';
    if ~isempty(Rb) && ~isempty(Rd)
    if strcmp(Config.Filter.RegMethod,'wnglimit')
        epsilon_f=testepsilons(1);
        w(:,ifreq)=getaccweights(Rd,Rb,I,epsilon_f,a_target);
        wngtest(1,ifreq)=getcurrWNG(w(:,ifreq),a_targets);
        wngbound=0;
        if wngtest(1,ifreq)<WNGConst(ifreq)-wngbound
        for etest=2:length(testepsilons)
            wtest=getaccweights(Rd,Rb,I,testepsilons(etest),a_target);
            wngtest(etest,ifreq)=getcurrWNG(wtest,a_targets);
        end
        [~,regind]=min(abs(WNGConst(ifreq)-wngtest(:,ifreq)));
        epsilon_f=testepsilons(regind);
        w(:,ifreq)=getaccweights(Rd,Rb,I,epsilon_f,a_target);
        end
        
    elseif strcmp(Config.Filter.RegMethod,'external')
        epsilon_f=epsilon(ifreq);
        w(:,ifreq)=getaccweights(Rd,Rb,I,epsilon_f,a_target);
    end
    end
    
end

ACCWeights.fd=w;

if strcmp(Config.Filter.FreqMode,'filter')
    
    ACCWeights.td=GetTDWeights(w,nMics,Config);
    
else
    ACCWeights.td=[];
end

end

function w=getaccweights(Rd,Rb,I,epsilon,a_target)

[w_hat,~]=eigs( (Rd+epsilon*I)\Rb ,1);

onaxisresp=w_hat'*a_target;

w=w_hat./onaxisresp;

%single target direction equates above expression to variation of SDB
%w=( (Rd+epsilon*I)\a_target) / ( (a_target'*((Rd+epsilon*I)\a_target)) );


end

