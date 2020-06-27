function [ArrayMan,Config]=GetArrayManifold(Config)
% function ArrayMan=GetArrayManifold(Config)
%
% Calculate the array manifold vector based on the parameters in Config.
%
% input argument: 
%   Config: configuration structure containing all the settings.

% output arguments: 
%   Config: updated configuration structure containing all the settings.
%   ArrayMan: struct containing the array manifold transfer functions.
%       
%       ArrayMan.Setup contains the ideal array manifold vector.
%       ArrayMan.Monitor contains the array manifold vector based on 
%           perturbed microphone positions (useful for robustness analysis)
%
%       For each Setup or Monitor a number of fields are output:
%       ArrayMan.a is the array manifold transfer function.
%       ArrayMan.Rvv contains the spatial correlation matrix between the
%           microphones.
%       ArrayMan.Rvv_ana is the analytical spatial correlation matrix.
%       ArrayMan.f is the frequency vector.
%       ArrayMan.LookAng.Az is the vector of steering angles in azimuth.
%       ArrayMan.LookAng.El is the vector of steering angles in elevation.



Config=checkConfigArrayMan(Config);
CurrConfig = Config;
% error catching
if ~strcmp(CurrConfig.ArrayMan.ManifoldType,'freefield') && ~strcmp(CurrConfig.ArrayMan.ManifoldType,'spherical_scatterer') ...
        && ~strcmp(CurrConfig.ArrayMan.ManifoldType,'cylindrical_scatterer')
    disp('choose from "freefield", "spherical_scatterer" or "cylindrical_scatterer" manifolds');
    error('incorrect manifold type');
end

CurrConfig.Array.CurrSelfNoise = 0;
CurrConfig.Array.CurrCalErrorMagdB = 0;
CurrConfig.Array.CurrCalErrorFreqdB = 0;
CurrConfig.ArrayMan.CurrDiffSound = 0;
CurrConfig.ArrayMan.CurrIsoNoise = 0;
ArrayMan.Setup=ArrayManifold(CurrConfig);

disp('setup manifold calculated');

poserrormag=CurrConfig.Array.Offset;
nnoise=length(poserrormag);
if nnoise
    calcnoise=1;
    OrigPos=CurrConfig.Array.MicPos;
else
    calcnoise=0;
end

if ~isfield(CurrConfig.ArrayMan,'Ntests'); CurrConfig.ArrayMan.Ntests = 1; end
ntests=CurrConfig.ArrayMan.Ntests;

if calcnoise
    ArrayMan.Monitor(nnoise+1,ntests)=struct('a',[],'Rvv',[],'Rvv_ana',[],'f',[],'LookAng',[],'OffsetDistance',[]);
    ArrayMan.Setup.OffsetDistance = []; 
    ArrayMan.Monitor(1,1)=ArrayMan.Setup;
    for inoise=1:nnoise
        for itest=1:ntests
        %spatially uncorrelated normally-distributed positioning error implemented here
        [OffsetPos, OffsetDistance] = ArrayOffset(OrigPos,CurrConfig.Array.Offset(inoise),CurrConfig.ArrayMan.ManifoldType,CurrConfig.Array.Geo);
        CurrConfig.Array.MicPos = OffsetPos;
        %diffuse sound magnitude
        CurrConfig.ArrayMan.CurrDiffSound=10.^(-CurrConfig.ArrayMan.DRRdB(inoise)/20);
        %isotropic noise magnitude
        CurrConfig.ArrayMan.CurrIsoNoise=10.^(-CurrConfig.ArrayMan.SNRIsodB(inoise)/20);
        %self noise magnitude from SNR value
        CurrConfig.Array.CurrSelfNoise=10^(-CurrConfig.Array.SNR(inoise)/20);
        %frequency independent magnitude error from dB error 
        CurrConfig.Array.CurrCalErrorMagdB=CurrConfig.Array.CalErrorMagdB(inoise);
        %magnitude dependent magnitude error from dB error
        CurrConfig.Array.CurrCalErrorFreqdB=CurrConfig.Array.CalErrorFreqdB(inoise);
        ArrayMan1 = ArrayManifold(CurrConfig);
        ArrayMan1.OffsetDistance = OffsetDistance;
        ArrayMan.Monitor(inoise+1,itest)=ArrayMan1; clear ArrayMan1
        fprintf('error manifold %d of %d calculated\n',inoise,nnoise);
        fprintf('test %d of %d\n',itest,ntests);
        end
    end
    ArrayMan.Monitor(1,1).poserrormag=0;
    for inoise=1:nnoise
        for itest=1:ntests
        ArrayMan.Monitor(inoise+1,itest).poserrormag=poserrormag(inoise);
        ArrayMan.Monitor(inoise+1,itest).SNR=CurrConfig.Array.SNR(inoise);
        ArrayMan.Monitor(inoise+1,itest).CalErrorMagdB=CurrConfig.Array.CalErrorMagdB(inoise);
        ArrayMan.Monitor(inoise+1,itest).CalErrorFreqdB=CurrConfig.Array.CalErrorFreqdB(inoise);
        ArrayMan.Monitor(inoise+1,itest).DRRdB=CurrConfig.ArrayMan.DRRdB(inoise);
        ArrayMan.Monitor(inoise+1,itest).SNRIsodB=CurrConfig.ArrayMan.SNRIsodB(inoise);
        end
    end
    
else
    ArrayMan.Monitor=ArrayMan.Setup;
    ArrayMan.Monitor.SNR = inf;
    ArrayMan.Monitor.CalErrorMagdB = 0;
    ArrayMan.Monitor.CalErrorFreqdB = 0;
    ArrayMan.Monitor.DRRdB = 0;
    ArrayMan.Monitor.SNRIsodB = inf;
    
end




end
function ArrayMan=ArrayManifold(Config)

LookAngEl = Config.ArrayMan.LookAng.El;
LookAngAz = Config.ArrayMan.LookAng.Az;
nlookangles=length(LookAngEl);

LookAngAzRad=deg2rad(LookAngAz);
LookAngElRad=deg2rad(LookAngEl);
% LookAngAz3D = reshape(repmat(LookAngAz,length(LookAngEl),1),nlookangles,1);
% LookAngEl3D = reshape(repmat(LookAngEl',1,length(LookAngAz)),nlookangles,1);
% LookAng3D = [LookAngAz3D LookAngEl3D];

nmics=size(Config.Array.MicPos,1);

if strcmp(Config.Filter.FreqMode,'discrete')
    f=Config.Filter.FreqVector;
else
    NfftNyquist = floor(Config.Filter.Nfft/2)+1;
    f = 0:Config.Filter.Fs/Config.Filter.Nfft:Config.Filter.Fs/Config.Filter.Nfft*(NfftNyquist-1);
end
nfreqs=length(f);

% array manifold (frequency domain)
% dimensions mics x look angles x frequency bins
a=zeros(nmics,nlookangles,nfreqs);

% spatial correlation matrix Rvv (frequency domain)
% dimensions mics x mics x frequency bins
R_vvangs=zeros(nmics,nmics,nlookangles); % tmp variable to store angles
R_vv=zeros(nmics,nmics,nfreqs); % Rvv
%R_vv_ana = zeros(nmics,nmics,length(f)); % Rvv analytical expresion for spherically/cylindrically isotropic sound fields

%calculating distances between each pair of microphones for analytical
%expression of noise coherence matrix
Mic_distances = zeros(nmics,nmics);
for n=1:nmics
    for m=1:nmics
        Mic_distances(n,m) = sqrt((Config.Array.MicPos(n,1)-Config.Array.MicPos(m,1))^2 + ...
            (Config.Array.MicPos(n,2)-Config.Array.MicPos(m,2))^2 + (Config.Array.MicPos(n,3)-Config.Array.MicPos(m,3))^2);
    end
end

k=2*pi*f./Config.ArrayMan.c;
r=Config.Array.MicPos;
r_norm2 = zeros(nmics,1);
for imics=1:nmics
    r_norm2(imics) = sqrt(r(imics,1)^2 + r(imics,2)^2 + r(imics,3)^2);
end
r_norm2 = mean(r_norm2);
rl_norm2 = Config.ArrayMan.PerformerDistance;
%coordinate system to match my simulations
[rlx,rly,rlz]=sph2car(LookAngAz,LookAngEl,Config.ArrayMan.PerformerDistance);
[rlux,rluy,rluz]=sph2car(LookAngAz,LookAngEl,1);
rl = cat(2,rlx,rly,rlz);
rlu = cat(2,rlux,rluy,rluz);

ErrCal = ones(nmics,1);
ErrResp = ones(nmics,1);
SelfNoise = zeros(nfreqs,nmics);
NoiseIso = zeros(nfreqs,nlookangles);
Diffsound = zeros(1,nlookangles);

if ~isempty(Config.ArrayMan.CurrDiffSound) && Config.ArrayMan.CurrDiffSound~=inf
    doDiff=1;
else
    doDiff=0;
end
if ~isempty(Config.ArrayMan.CurrIsoNoise) && Config.ArrayMan.CurrIsoNoise~=inf
    doIso=1;
else
    doIso=0;
end
%if any of the errors are to be modelled
if ~isempty(Config.Array.Offset)
    %frequency independent calibration error
        ErrCal=10.^(-Config.Array.CurrCalErrorMagdB.*randn([nmics,1])/20);
        %sensor self noise as spatially uncorrelated as a function of
        %microphone
        
        for imic=1:nmics
            SelfNoise(:,imic) = Config.Array.CurrSelfNoise/sqrt(2).*(randn(nfreqs,1)+1i*randn(nfreqs,1));
        end
        
        %isotropic noise as spatially uncorrelated as a function of
        %angle but correlated as a function of microphone. 
        if doDiff; Diffsound=Config.ArrayMan.CurrDiffSound.*exp(1i*2*pi*rand(1,nlookangles)); end

        %frequency dependent calibration error
        ErrResp=10.^(-Config.Array.CurrCalErrorFreqdB.*randn([nmics,1])/20);
end

% calculate array manifold based on free field or spherically diffracted distances
if strcmp(Config.ArrayMan.ManifoldType,'freefield')
        
    if strcmp(Config.ArrayMan.LookDistance,'focalpoint')
        % for focal point, calculate the positions of focal points in the look directions
        
        % calculate the distances from each microphone to a source at the
        % focal point and far field source
        R=zeros(nmics,nlookangles);
        Ru=zeros(nmics,nlookangles);

        for ilook=1:nlookangles
            R(:,ilook)=sqrt(abs(rlx(ilook)-r(:,1)).^2 ...
                + abs(rly(ilook)-r(:,2)).^2  ...
                + abs(rlz(ilook)-r(:,3)).^2);
            Ru(:,ilook)=sqrt(abs(rlux(ilook)-r(:,1)).^2 ...
                + abs(rluy(ilook)-r(:,2)).^2  ...
                + abs(rluz(ilook)-r(:,3)).^2);
        
            if doIso; NoiseIso(:,ilook) = Config.ArrayMan.CurrIsoNoise/sqrt(2).*(randn(nfreqs,1)+1i*randn(nfreqs,1)); end
            
        end
        for ifreq=1:nfreqs  
            SelfNoisef=SelfNoise(ifreq,:);     
            for ilook=1:nlookangles
                % Array manifold vector steered in each direction
                Rl=R(:,ilook);   
                pdir = 1./(4.*pi.*Rl).*exp(1i*k(ifreq)*Rl);
                %isotropic noise
                Rul=Ru(:,ilook);   
                pdiff = exp(1i*k(ifreq)*Rul);
                atmp= ErrCal.*ErrResp.*(pdir + 1/nlookangles.*((Diffsound(ilook) + NoiseIso(ifreq,ilook)).*pdiff)) + SelfNoisef.';

                a(:,ilook,ifreq)=atmp;
                % isotropic acoustical noise
                R_vvangs(:,:,ilook)=atmp*atmp'; % (this calculation MUCH quicker with atmp)
            end
            R_vv(:,:,ifreq)=sum(R_vvangs.*shiftdim(Config.ArrayMan.Qweights,-2),3);
        end
        
        if strcmp(Config.ArrayMan.SteerSpace,'2Daz') || strcmp(Config.ArrayMan.SteerSpace,'2Del')
            %for Matlab 2017
            %R_vv_ana = besselj(0, 2*pi.*Mic_distances.*shiftdim(f,-1)./Config.ArrayMan.c);
            %for older Matlabs dimensions need to be the same
            R_vv_ana = besselj(0, 2*pi.*repmat(Mic_distances,1,1,nfreqs).*repmat(shiftdim(f,-1),nmics,nmics,1)./Config.ArrayMan.c);
        else
            %for Matlab 2017
            %R_vv_ana = sinc(2.*Mic_distances.*shiftdim(f,-1)./Config.ArrayMan.c); % pi already included in sinc function
            %for older Matlabs dimensions need to be the same
            R_vv_ana = besselj(0, 2.*repmat(Mic_distances,1,1,nfreqs).*repmat(shiftdim(f,-1),nmics,nmics,1)./Config.ArrayMan.c);
        end
        
    elseif strcmp(Config.ArrayMan.LookDistance,'farfield')
       
         if doIso   
             for ilook=1:nlookangles
                NoiseIso(:,ilook) = Config.ArrayMan.CurrIsoNoise/sqrt(2).*(randn(nfreqs,1)+1i*randn(nfreqs,1));
             end
         end
         for ifreq=1:nfreqs
             SelfNoisef=SelfNoise(ifreq,:);
            for ilook=1:nlookangles
                % Array manifold vector steered in each direction
                pdir = exp(1i*k(ifreq)*r*shiftdim(rlu(ilook,:),1));
                
                atmp= ErrCal.*ErrResp.*pdir.*(1 + 1/nlookangles.*((Diffsound(ilook) + NoiseIso(ifreq,ilook)))) + SelfNoisef.';                   
                a(:,ilook,ifreq)=atmp;
                % cross- and auto- power spectral densities
%                     cross = atmp*atmp';
%                     auto = diag(cross);
                R_vvangs(:,:,ilook)=atmp*atmp'; % (this calculation MUCH quicker with atmp)
                %R_vvangs(:,:,ilook,ifreq)=cross./sqrt(auto*auto'); % (this calculation MUCH quicker with atmp)

            end                
            R_vv(:,:,ifreq)=sum(R_vvangs.*shiftdim(Config.ArrayMan.Qweights,-2),3);
        end

        if strcmp(Config.ArrayMan.SteerSpace,'2Daz') || strcmp(Config.ArrayMan.SteerSpace,'2Del')
            %for Matlab 2017
            %R_vv_ana = besselj(0, 2*pi.*Mic_distances.*shiftdim(f,-1)./Config.ArrayMan.c);
            %for older Matlabs dimensions need to be the same
            R_vv_ana = besselj(0, 2*pi.*repmat(Mic_distances,1,1,nfreqs).*repmat(shiftdim(f,-1),nmics,nmics,1)./Config.ArrayMan.c);
        else
            %for Matlab 2017
            %R_vv_ana = sinc(2.*Mic_distances.*shiftdim(f,-1)./Config.ArrayMan.c); % pi already included in sinc function
            %for older Matlabs dimensions need to be the same
            R_vv_ana = besselj(0, 2.*repmat(Mic_distances,1,1,nfreqs).*repmat(shiftdim(f,-1),nmics,nmics,1)./Config.ArrayMan.c);
        end
        
    else
        error('incorrect manifold look distance');
    end
    
elseif strcmp(Config.ArrayMan.ManifoldType,'spherical_scatterer')
    R_vv_ana = [];
    for ilook=1:nlookangles
        if strcmp(Config.ArrayMan.LookDistance,'focalpoint') 
            GN = Neumann_Greens_function_v2(r,rl(ilook,:),r_norm2,rl_norm2,k);
            a(:,ilook,:)=reshape(GN,size(GN,1),1,size(GN,2));
            if doDiff
                pdiff = pressure_on_rigid_sphere_from_plane_wave(r,rlu(ilook,:),k);
                adiff(:,ilook,:)=reshape(pdiff,size(pdiff,1),1,size(pdiff,2));
            end
        elseif strcmp(Config.ArrayMan.LookDistance,'farfield')  
            pdiff = pressure_on_rigid_sphere_from_plane_wave(r,rlu(ilook,:),k);
            a(:,ilook,:) = reshape(pdiff,size(pdiff,1),1,size(pdiff,2));
            adiff(:,ilook,:) = a(:,ilook,:);
        else  
            error('incorrect manifold look distance');
        end 
        if doIso; NoiseIso(:,ilook) = Config.ArrayMan.CurrIsoNoise/sqrt(2).*(randn(nfreqs,1)+1i*randn(nfreqs,1)); end
    end
        
    for ifreq=1:nfreqs
        SelfNoisef=SelfNoise(ifreq,:);
        for ilook=1:nlookangles
            %isotropic noise
            atmp = ErrCal.*ErrResp.*a(:,ilook,ifreq) + adiff(:,ilook,ifreq).*(1/nlookangles.*(Diffsound(ilook) + NoiseIso(ifreq,ilook))) + SelfNoisef.';
            a(:,ilook,ifreq)=atmp;
            R_vvangs(:,:,ilook)=atmp*atmp';
        end
        R_vv(:,:,ifreq)=sum(R_vvangs.*shiftdim(Config.ArrayMan.Qweights,-2),3);
    end        
   
elseif strcmp(Config.ArrayMan.ManifoldType,'cylindrical_scatterer')
    R_vv_ana = [];
    for ilook=1:nlookangles
        if strcmp(Config.ArrayMan.LookDistance,'focalpoint') 
            error('point source array manifold for cylindrical scatter not implemented yet');
        elseif strcmp(Config.ArrayMan.LookDistance,'farfield') 
            pdiff = pressure_on_rigid_cylinder_from_plane_wave(r,rlu(ilook,:),k);
            a(:,ilook,:) = reshape(pdiff,size(pdiff,1),1,size(pdiff,2));
            adiff(:,ilook,:) = a(:,ilook,:);
        else  
            error('incorrect manifold look distance');
        end 
        if doIso; NoiseIso(:,ilook) = Config.ArrayMan.CurrIsoNoise/sqrt(2).*(randn(nfreqs,1)+1i*randn(nfreqs,1)); end
    end
    for ifreq=1:nfreqs
        SelfNoisef=SelfNoise(ifreq,:);
        for ilook=1:nlookangles
            %isotropic noise
            atmp = ErrCal.*ErrResp.*a(:,ilook,ifreq) + adiff(:,ilook,ifreq).*(1/nlookangles.*(Diffsound(ilook) + NoiseIso(ifreq,ilook))) + SelfNoisef.';
            a(:,ilook,ifreq)=atmp;
            R_vvangs(:,:,ilook)=atmp*atmp';
        end
        R_vv(:,:,ifreq)=sum(R_vvangs.*shiftdim(Config.ArrayMan.Qweights,-2),3);
    end        
    
end 

    ArrayMan.a=a;
    ArrayMan.Rvv=R_vv;
    ArrayMan.Rvv_ana = R_vv_ana;
    ArrayMan.f=f;
    ArrayMan.LookAng.Az=LookAngAz;
    ArrayMan.LookAng.El=LookAngEl;

end
