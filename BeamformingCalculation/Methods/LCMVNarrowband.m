function [ LCMVWeights ] = LCMVNarrowband(Config,ArrayMan,R_xx)
% function [ LCMVWeights ] = LCMVNarrowband(Config,ArrayMan,R_xx)
% Calculate narrowband filter weights for LCMV (generalized MVDR), based
% on the parameters specified in the Config struct.
%
% input arguments:
%   Config: configuration struct containing all the settings.
%   ArrayMan: struct containing the array manifold transfer functions.
%
% output arguments:
%   LCMVWeights: beamforming weights struct. 

nFreqs=length(ArrayMan.f);
nMics=length(Config.Array.MicPos);
w=zeros(nMics,nFreqs);
I=eye(nMics);

targetAz = Config.Filter.CurrTarget(2);
targetEl = Config.Filter.CurrTarget(1);
intAngles = Config.Filter.CurrInts;
nints=length(intAngles)/2;
[~,~,targetInd] = getNearestAngle(Config,targetEl,targetAz);
intInd = zeros(nints,1);
for iint=1:nints
    [~,~,intInd(iint)] = getNearestAngle(Config,intAngles(iint*2-1),intAngles(iint*2));
end

% constraint vector
g=db2mag(-40)*ones(1+nints,1);
g(1)=1; % assume a single target location

testepsilons=logspace(-6,3,200);         
wngtest=zeros(length(testepsilons),nFreqs);
for ifreq=1:nFreqs
    
    % array manifold in target and notch locations
    
    %C_k=squeeze(ArrayMan.a(:,[targetElind,intElinds],[targetAzind,intAzinds],ifreq));
    %C_k=reshape(C_k,size(C_k,1),size(C_k,2)*size(C_k,3));
    
    Ctar = squeeze(ArrayMan.a(:,targetInd,ifreq));
    if ~isempty(intInd)
        Cint = squeeze(ArrayMan.a(:,intInd,ifreq));
        C_k=[Ctar,Cint];
    else
        C_k = Ctar;
    end
    
    R_xxk=(R_xx(:,:,ifreq));
    if sum(sum(~isnan(C_k)))==numel(C_k)
    if strcmp(Config.Filter.RegMethod,'wnglimit')
        %epsilon_f=0.0001;
        epsilon_f=testepsilons(1);
        w(:,ifreq)=getlcmvweights(R_xxk,C_k,g,I,epsilon_f);
        wngtest(1,ifreq)=getcurrWNG(w(:,ifreq),Ctar);
        wngbound=0;
        if wngtest(1,ifreq)<Config.Filter.RegParam-wngbound
        for etest=2:length(testepsilons)
            wtest=getlcmvweights(R_xxk,C_k,g,I,testepsilons(etest));
            wngtest(etest,ifreq)=getcurrWNG(wtest,Ctar);
        end
        [~,regind]=min(abs(Config.Filter.RegParam-wngtest(:,ifreq)));
        epsilon_f=testepsilons(regind);
        w(:,ifreq)=getlcmvweights(R_xxk,C_k,g,I,epsilon_f);
        end
        
    elseif strcmp(Config.Filter.RegMethod,'external')
        epsilon_f=Config.Filter.RegParam;
        w(:,ifreq)=getlcmvweights(R_xxk,C_k,g,I,epsilon_f);
    end
    LCMVWeights.epsilon(ifreq)=epsilon_f;
    end
    
end

LCMVWeights.fd=w;

if strcmp(Config.Filter.FreqMode,'filter')
    
    LCMVWeights.td=GetTDWeights(w,nMics,Config);
        
else
    LCMVWeights.td=[];
end

end

function w=getlcmvweights(R_xxk,C_k,g,I,epsilon_f)

w=((R_xxk+epsilon_f*I)\C_k) * ( (C_k'*((R_xxk+epsilon_f*I)\C_k)) \ g);
end
