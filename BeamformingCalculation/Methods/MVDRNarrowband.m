function [ MVDRWeights ] = MVDRNarrowband(Config,Audio)
% function [ MVDRWeights ] = MVDRNarrowband(Config,Audio)
% Calculate narrowband filter weights for MVDR, based on the parameters 
% specified in the Config struct.
%
% input arguments:
%   Config: configuration struct containing all the settings.
%   ArrayMan: struct containing the array manifold transfer functions.
%
% output arguments:
%   MVDRWeights: beamforming weights struct. 


% distances between the far field source and the microphone positions
[Srcx,Srcy]=pol2cart(deg2rad(Config.Angle),Config.FarFieldDist);

SegSamples=1000;
DataSize=size(Audio);
StartSamples=1:SegSamples:max(DataSize);

R_x=zeros(length(StartSamples),Config.Nfft,min(DataSize),min(DataSize));
%tic
for j=2:length(StartSamples)
    AudioSeg=Audio(:,StartSamples(j-1):StartSamples(j)-1);
    DataFD=fft(AudioSeg,Config.Nfft,2);
    
    for k=1:Config.Nfft/2
        R_x(j,k,:,:)=DataFD(:,k)*DataFD(:,k)';
    end
end
%toc
E_R_x=squeeze(mean(R_x));
f=Config.f;

Delayx=zeros(length(Srcx),length(Config.MicPos));
Delayy=zeros(length(Srcy),length(Config.MicPos));
Delay=zeros(size(Delayx));
w=zeros(Config.Nfft,length(Config.MicPos),length(Srcx));
d=zeros(size(w));
for i=1:length(Srcx)
    Delayx(i,:)=((1/Config.c)*(abs(Srcx(i)-Config.MicPos(:,1))));
    Delayy(i,:)=((1/Config.c)*(abs(Srcy(i)-Config.MicPos(:,2))));
    Delay(i,:)=sqrt(Delayx(i,:).^2 + Delayy(i,:).^2);
    
    for k=2:Config.Nfft/2
        d(k,:,i)=exp(-1i*2*pi*f(k)*Delay(i,:));
        w(k,:,i)=(squeeze(E_R_x(k,:,:))\d(k,:,i)')/(d(k,:,i)*(squeeze(E_R_x(k,:,:))\d(k,:,i)'));
    end
    w(Config.Nfft/2 +2:end,:,i)=flipud(conj(w(2:Config.Nfft/2,:,i)));
    
end
MVDRWeights=real(ifft(w));

end

