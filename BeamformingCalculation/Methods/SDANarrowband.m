function [ SDAWeights ] = SDANarrowband(Config,ArrayMan)
% function [ SDAWeights ] = SDANarrowband(Config,ArrayMan)
%
% Calculate narrowband filter weights for the superdirective array (SDA), 
% based on the parameters specified in the Config struct.
% (Bai et al., Acoustic Array Systems: theory, implementation and
% application, Wiley, 2013, p.122)
%
% input arguments: 
%   Config: configuration struct containing all the settings.
%   ArrayMan: struct containing the array manifold transfer functions.
%
% output arguments: 
%   SDAWeights: beamforming weights struct. 


nFreqs=length(ArrayMan.f);
nMics=length(Config.Array.MicPos);
w=zeros(nMics,nFreqs);
I=eye(nMics);

targetAz = Config.Filter.CurrTarget(2);
targetEl = Config.Filter.CurrTarget(1);
[~,~,targetInd] = getNearestAngle(Config,targetEl,targetAz);

testepsilons=logspace(-6,3,200);         
wngtest=zeros(length(testepsilons),nFreqs);
for ifreq=1:nFreqs
    
    % Array manifold for that frequency and direction
    a_target=ArrayMan.a(:,targetInd,ifreq);
    Rvv=ArrayMan.Rvv(:,:,ifreq);
    %epsilon(ifreq)=1;%eigs(Rvv,1)./Config.Filter.RegParam;
    if sum(~isnan(a_target),1)==size(a_target,1) || sum(sum(~isnan(Rvv)))==numel(Rvv)
    if strcmp(Config.Filter.RegMethod,'wnglimit')
        %epsilon_f=0.0001;
        epsilon_f = testepsilons(1);
        w(:,ifreq)=getsdaweights(a_target,Rvv,I,epsilon_f);
        wngtest(1,ifreq)=getcurrWNG(w(:,ifreq),a_target);
        wngbound=0;
        if wngtest(1,ifreq)<Config.Filter.RegParam-wngbound
        for etest=2:length(testepsilons)
            wtest=getsdaweights(a_target,Rvv,I,testepsilons(etest));
            wngtest(etest,ifreq)=getcurrWNG(wtest,a_target);
        end
        [~,regind]=min(abs(Config.Filter.RegParam-wngtest(:,ifreq)));
        epsilon_f=testepsilons(regind);
        w(:,ifreq)=getsdaweights(a_target,Rvv,I,epsilon_f);
        end
        
    elseif strcmp(Config.Filter.RegMethod,'external')
        epsilon_f=Config.Filter.RegParam;
        w(:,ifreq)=getsdaweights(a_target,Rvv,I,epsilon_f);
    end
    SDAWeights.epsilon(ifreq)=epsilon_f;
    end
    
    
    
    
    
end

SDAWeights.fd=w;

if strcmp(Config.Filter.FreqMode,'filter')
    
    SDAWeights.td=GetTDWeights(w,nMics,Config);
        
else
    SDAWeights.td=[];
end

end

function w=getsdaweights(a_target,Rvv,I,epsilon_f)

w=( (Rvv+epsilon_f*I)\a_target) / ( (a_target'*((Rvv+epsilon_f*I)\a_target)) );
    
end
