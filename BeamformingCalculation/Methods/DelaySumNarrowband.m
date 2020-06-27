function [ DSWeights ] = DelaySumNarrowband(Config,ArrayMan)
% function [ DSWeights ] = DelaySumNarrowband(Config,ArrayMan)
% Calculate narrowband filter weights for delay and sum beamforming, based
% on the parameters specified in the Config struct.
%
% input arguments:
%   Config: configuration struct containing all the settings.
%   ArrayMan: struct containing the array manifold transfer functions.
%
% output arguments:
%   DSWeights: beamforming weights struct. 

nFreqs=length(ArrayMan.f);
nMics=length(Config.Array.MicPos);
w=zeros(nMics,nFreqs);

targetAz = Config.Filter.CurrTarget(2);
targetEl = Config.Filter.CurrTarget(1);
[~,~,targetInd] = getNearestAngle(Config,targetEl,targetAz);

for k=1:nFreqs
    w(:,k)=1/nMics.*ArrayMan.a(:,targetInd,k);
end


DSWeights.fd=w;

if strcmp(Config.Filter.FreqMode,'filter')
    
    DSWeights.td=GetTDWeights(w,nMics,Config);
        
else
    DSWeights.td=[];
end


end





