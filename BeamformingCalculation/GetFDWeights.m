function fdweights=GetFDWeights(w,Config)
% fdweights=GetFDWeights(w,Config)
% function to compute the beamforming weights in frequency domain
% from the time domain counterpart.
%
% input arguments:
%   w: beamforming weights in the time domain.
%   Config: configuration structure containing all the settings.
%
% output arguments:
%   fdweights: beamforming weights in the frequency domain. 

nfft=Config.Filter.Nfft;
filt=fft(ifftshift(w,1),[],1);
%complex conjugation transposition to get w rather than w^H
fdweights = filt(1:floor(nfft/2+1),:)';

end