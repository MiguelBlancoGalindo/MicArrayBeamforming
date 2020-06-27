function [ MNWeights ] = MNNarrowband(Config,ArrayMan)
% function [ MNWeights ] = MNNarrowband(Config,ArrayMan)
% Calculate narrowband filter weights for minimum norm beamformer, based
% on the parameters specified in the Config struct.
%
% input arguments:
%   Config: configuration struct containing all the settings.
%   ArrayMan: struct containing the array manifold transfer functions.
%
% output arguments:
%   MNWeights: beamforming weights struct. 

nFreqs=length(ArrayMan.f);
nMics=length(Config.Array.MicPos);
w=zeros(nMics,nFreqs);

phi_t = Config.Filter.CurrTarget(2);
theta_t = Config.Filter.CurrTarget(1);
[~,~,itar] = getNearestAngle(Config,theta_t,phi_t);
    
d_d = getTargetPattern(Config);

% response vector from nulls and target response
inulls = findNulls(d_d);
nnulls = length(inulls); 
d = [1 zeros(1,nnulls)];
I=eye(nnulls+1);

testepsilons=logspace(-6,3,200);   
wngtest=zeros(length(testepsilons),nFreqs);
d_look=zeros(nFreqs,1);

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

for ifreq=1:nFreqs
    
    % array manifold matrix
    A_int=ArrayMan.a(:,inulls,ifreq);
    a_tar=ArrayMan.a(:,itar,ifreq);
    A=cat(2,a_tar,A_int);
    % check for potential nan in array manifold
    nanPos = find(sum(isnan(A),1)==size(A,1));
    dfreq = d;
    dfreq(nanPos) = [];
    % removing potential nan direction
    A(:,nanPos) = [];
    epsilon_f=testepsilons(1);
    if strcmp(Config.Filter.RegMethod,'wnglimit') && ~isempty(A)

        w(:,ifreq)=getmnweights(A,dfreq,I,epsilon_f);
        wngtest(1,ifreq)=getcurrWNG(w(:,ifreq),a_tar);
        wngbound=0;
        if wngtest(1,ifreq)<WNGConst(ifreq)-wngbound
        for etest=2:length(testepsilons)
            wtest=getmnweights(A,dfreq,I,testepsilons(etest));
            wngtest(etest,ifreq)=getcurrWNG(wtest,a_tar);
        end
        [~,regind]=min(abs(WNGConst(ifreq)-wngtest(:,ifreq)));
        epsilon_f=testepsilons(regind);
        w(:,ifreq)=getmnweights(A,dfreq,I,epsilon_f);
        
        end
        
    elseif strcmp(Config.Filter.RegMethod,'external')
        
        epsilon_f=epsilon(ifreq);
        w(:,ifreq)=getmnweights(A,dfreq,I,epsilon_f);

    end
    d_look(ifreq) = w(:,ifreq)'*a_tar;
    MNWeights.epsilon(ifreq)=epsilon_f;
    
    
end

%equalising response at look direction
%finding last frequency whose response is less than 0.01 dB from flat
%region
%first calculate the positions of d_look whose gradient is less than thresh
threshdB = 0.01;
ipos = find(abs(diff(db(d_look)))<threshdB);
%then get gradient of those positions to find out where the largest flat
%region is
grad = diff(ipos);
%invert the gradient for non consecutive values before taking cumulative
%sum (from beginning to end)
grad(grad>1) = -grad(grad>1);
[~,imax] = max(cumsum(grad));
ifo = ipos(imax)+2; %adding 2 because two gradient functions were used
h_eq = ones(nFreqs,1);
h_eq(1:ifo) = 1./abs(d_look(1:ifo));
w = w.*h_eq.';


MNWeights.fd=w;

if strcmp(Config.Filter.FreqMode,'filter')
    
    MNWeights.td=GetTDWeights(w,nMics,Config);
    
else
    MNWeights.td=[];
end

end

function w=getmnweights(A,d,I,epsilon_f)

w=A/(A'*A + epsilon_f*I)*d';

end

function iint = findNulls(d_d)

d_d = -abs(d_d);

[~,iint] = findpeaks(d_d); 

end

