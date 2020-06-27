function atd=GetTDArrayMan(a,nMics,Config)
% function atd=GetTDArrayMan(a,nMics,Config)
% function to compute the array manifold transfer functions in time domain
% from the frequency domain counterpart.
%
% input arguments:
%   a: array manifold transfer function in the frequency domain.
%   nMics: number of microphones.
%   Config: configuration structure containing all the settings.
%
% output arguments:
%   atd: array manifold transfer function in the time domain. 

nfft=Config.Filter.Nfft;
filt=zeros(nMics,nfft);

filt(:,1:nfft/2+1)=a;
filt(:,1)=zeros(nMics,1);
filt(:,nfft/2+1)=zeros(nMics,1);
filt(:,nfft/2+2:end)=fliplr(conj(a(:,2:nfft/2)));

%Transposition applied to match dimensions of signals later on
atd=real(ifftshift(ifft(filt.',[],1),1));

%Apply 50% Tukey window to ensure decay at the beginning and end
atd = atd.*tukeywin(nfft,0.5);