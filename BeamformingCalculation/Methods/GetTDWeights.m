function tdweights=GetTDWeights(w,nMics,Config)
% function tdweights=GetTDWeights(w,nMics,Config)
% function that calculates the time-domain weights from the equivalent
% frequency-domain counterparts. 
%
% input arguments: 
%   w: microphone array filter weights in the frequency domain.
%   nMics: number of microphones
%   Config: configuration struct containing all the settings.
%
% output arguments: 
%   tdweights: time-domain weights. 

nfft=Config.Filter.Nfft;
filt=zeros(nMics,nfft);

filt(:,1:nfft/2+1)=w;
filt(:,1)=zeros(nMics,1);
filt(:,nfft/2+1)=zeros(nMics,1);
filt(:,nfft/2+2:end)=fliplr(conj(w(:,2:nfft/2)));

%Weights are conjugate-transposed before taking the inverse Fourier 
%Transform. 
tdweights=real(ifftshift(ifft(filt',[],1),1));
%Apply 50% Tukey window to ensure decay at the beginning and end
tdweights = tdweights.*tukeywin(nfft,0.5);
end