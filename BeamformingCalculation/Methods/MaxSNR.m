function [ MaxSNRWeights ] = MaxSNR(Config)
% Calculate narrowband filter weights for delay and sum beamforming, based
% on the parameters specified in the Config struct

% distances between the far field source and the microphone positions
[Srcx,Srcy]=pol2cart(deg2rad(1:360),Config.FarFieldDist);

TargetAngleInd=zeros(1,length(Config.Angle));
for j=1:length(Config.Angle)
    if Config.Angle(j)<1
        TargetAngleInd(j)=360+Config.Angle;
    else
        TargetAngleInd(j)=Config.Angle;
    end
end
passBeam=10;
stopBeam=20;
%R_vvangs=zeros(length(Config.MicPos),length(Config.MicPos),length(Srcx));
f=((0:Config.Nfft-1)./Config.Nfft)*Config.Fs;
Delayx=zeros(length(Srcx),length(Config.MicPos));
Delayy=zeros(length(Srcy),length(Config.MicPos));
Delay=zeros(size(Delayx));
w=zeros(Config.Nfft,length(Config.MicPos),length(Config.Angle));
d=zeros(length(Config.MicPos),length(Srcx));
I=eye(length(Config.MicPos),length(Config.MicPos));
beta=100;

for i=1:360
    Delayx(i,:)=((1/Config.c)*(abs(Srcx(i)-Config.MicPos(:,1))));
    Delayy(i,:)=((1/Config.c)*(abs(Srcy(i)-Config.MicPos(:,2))));
    Delay(i,:)=sqrt(Delayx(i,:).^2 + Delayy(i,:).^2);
end

for j=1:length(Config.Angle)
    for k=2:Config.Nfft/2
        for i=1:360
            % Array manifold vector steered in each direction
            d(:,i)=exp(1i*2*pi*f(k)*Delay(i,:));
            %R_vvangs(:,:,i)=(d(:,i))*(d(:,i))';
        end
        % Array manifold vector steered at the target
        %passband=TargetAngleInd(j);
        passband = mod(((TargetAngleInd(j)-passBeam):(TargetAngleInd(j)+passBeam))-1, 360)+1;
        %R_s=sum(R_vvangs(:,:,passband),3)./size(R_vvangs(:,:,passband),3);
        d_target=d(:,passband);
        R_s=d_target*d_target';
        
        notstopInd = mod(((TargetAngleInd(j)-stopBeam):(TargetAngleInd(j)+stopBeam))-1, 360)+1;
        isStop = ones(360,1);
        isStop(notstopInd) = 0;
        stopband = isStop == 1;
        %R_n=sum(R_vvangs(:,:,stopband),3)./size(R_vvangs(:,:,stopband),3);
        R_n=d(:,stopband)*d(:,stopband)';
        W=(R_n+beta*I)\(R_s);
        [maxw,~]=eigs(W,1);
        w(k,:,j)=maxw/sqrt(length(W));
        
        
    end
    w(Config.Nfft/2 +2:end,:,j)=flipud(conj(w(2:Config.Nfft/2,:,j)));
    
end
MaxSNRWeights=real((ifft(w)));
%figure;plot(MaxSNRWeights);
end

