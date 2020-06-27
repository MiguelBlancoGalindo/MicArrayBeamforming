function [ Ainverse ] = SteeringMatrixPseudoInverse(Config,ArrayMan)
% function [ Ainverse ] = SteeringMatrixPseudoInverse(Config,ArrayMan)
% Calculate the pseudoinverse of the steering matrix, based
% on the parameters specified in the Config struct,(e.g. regularisation 
% parameter).
%
% input arguments:
%   Config: configuration struct containing all the settings.
%   ArrayMan: struct containing the array manifold steering matrix.
%
% output arguments:
%   Ainverse: pseudoinverse of the steering matrix. 

nFreqs=length(ArrayMan.f);
nMics=length(Config.Array.MicPos);
A_dagger=zeros(size(ArrayMan.a));
I=eye(nMics);    

testepsilons=logspace(-6,3,200);   
norm_Adagger_test=zeros(length(testepsilons),nFreqs);


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
    A=ArrayMan.a(:,:,ifreq);
    % check for potential nan in array manifold
    %nanPos = find(sum(isnan(A),1)==size(A,1));
    % removing potential nan direction
    %A(:,nanPos) = [];
    if strcmp(Config.Filter.RegMethod,'wnglimit')
        epsilon_f=testepsilons(1);
        Af_dagger=getpseudoinverse(A,I,epsilon_f);
        norm_Adagger_test(1,ifreq) = getcurrAnorm(A,Af_dagger);
        wngbound=0;
        if norm_Adagger_test(1,ifreq)>-(WNGConst(ifreq)-wngbound)
        for etest=2:length(testepsilons)
            Af_dagger_test=getpseudoinverse(A,I,epsilon_f);
            norm_Adagger_test(etest,ifreq) = getcurrAnorm(A,Af_dagger_test);
        end
        [~,regind]=min(abs(WNGConst(ifreq)-norm_Adagger_test(:,ifreq)));
        epsilon_f=testepsilons(regind);
        A_dagger(:,:,ifreq)=getpseudoinverse(A,I,epsilon_f);
        
        end
        
    elseif strcmp(Config.Filter.RegMethod,'external')
        epsilon_f=epsilon(ifreq);
        A_dagger(:,:,ifreq)=getpseudoinverse(A,I,epsilon_f);

    end
    Ainverse.epsilon(ifreq)=epsilon_f;
    
end

Ainverse.A_dagger=A_dagger;


end

function Adagger=getpseudoinverse(A,I,epsilon_f)

Adagger=(A*A' + epsilon_f*I)\A;

end
 
function norm_Adagger = getcurrAnorm(A,A_dagger)

norm_Adagger = db(norm(A_dagger)./norm(A));

end

