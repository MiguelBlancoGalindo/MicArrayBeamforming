function [Eval]=BFWeightEvaluation(BFWeights,Config,ArrayMan)
% function [Eval]=BFWeightEvaluation(BFWeights,Config,ArrayMan)
% calculate a range of evaluation metrics from the beamforming calculation
%
% input arguments:
%   BFWeights: beamforming weights struct
%   Config: configuration struct containing all the settings.
%   ArrayMan: struct containing the array manifold transfer functions.
% output arguments: 
%   Eval: evaluation struct containing several performance metrics. 


azTarget = Config.Filter.CurrTarget(2);
elTarget = Config.Filter.CurrTarget(1);

if ~isempty(Config.Filter.CurrInts)
    doInterf = true;
    azInt = Config.Filter.CurrInts(2:2:end);
    elInt = Config.Filter.CurrInts(1:2:end);
else
    doInterf = false;
    azInt = []; elInt = [];
end
nints=size(Config.Filter.CurrInts,2)/2;

[az,iaz] = sort(ArrayMan(1).LookAng.Az(ArrayMan(1).LookAng.El==0));
el1 = ArrayMan(1).LookAng.El(ArrayMan(1).LookAng.Az==0);
el2 = -ArrayMan(1).LookAng.El(ArrayMan(1).LookAng.Az==180) + 180;
el2(el2>180) = el2(el2>180)-360;
[el,iel] = sort([el1; el2]);

naz = length(az);
nel = length(el);
[~,iazTarget] = min(abs(az-azTarget));
[~,ielTarget] = min(abs(el-elTarget));
if doInterf
    [~,iazInt] = min(abs(az-azInt));
    [~,ielInt] = min(abs(el-elInt));
end
iazDirTarget = iazTarget;
ielDirTarget = ielTarget;
[~,~,iTarget] = getNearestAngle(Config,elTarget,azTarget);
iInt = zeros(nints,1);
for ii=1:nints
    [~,~,iInt(ii)] = getNearestAngle(Config,elInt(ii),azInt(ii));
end

[~,iaz0] = min(abs(az-0));
[~,iaz270] = min(abs(az-(-90)));
[~,iaz90] = min(abs(az-90));
iel0 = iaz0;
[~,ielminus90] = min(abs(el-(-90)));
[~,iel90] = min(abs(el-90));

iDirTarget = iTarget;
if ~isfield(Config.Scene,'AngRange') 
    azRange=60; elRange=60; 
else
    azRange=Config.Scene.AngRange(2); 
    elRange=Config.Scene.AngRange(1);
end

iRangeOn = find(azTarget-azRange/2<=ArrayMan(1).LookAng.Az & azTarget+azRange/2>=ArrayMan(1).LookAng.Az & ...
            elTarget-elRange/2<=ArrayMan(1).LookAng.El & elTarget+elRange/2>=ArrayMan(1).LookAng.El);
iRangeOff = find((azTarget-azRange/2>ArrayMan(1).LookAng.Az | azTarget+azRange/2<ArrayMan(1).LookAng.Az) & ...
    (elTarget-elRange/2>ArrayMan(1).LookAng.El | elTarget+elRange/2<ArrayMan(1).LookAng.El));
[~,~,iazRangemin] = getNearestAngle(Config,0,azTarget-azRange/2);
[~,~,iazRangemax] = getNearestAngle(Config,0,azTarget+azRange/2);
[~,~,ielRangemin] = getNearestAngle(Config,elTarget-elRange/2,0);
[~,~,ielRangemax] = getNearestAngle(Config,elTarget-elRange/2,0);

nfreqs=length(ArrayMan(1).f);
nlook = length(ArrayMan(1).LookAng.Az);

if strcmp(Config.ArrayMan.SteerSpace,'2Daz')
    Respdiraz=zeros(nfreqs,naz);
    Respshiftaz=zeros(nfreqs,naz);
    Respaz=zeros(nfreqs,naz);
    RespOffTarAz = zeros(nfreqs,1);
    DIaz=zeros(nfreqs,1);
    Contrastaz=zeros(size(DIaz));
    ContrastRangeaz=zeros(size(DIaz));
    WNG=zeros(size(DIaz));
    SLaz=zeros(size(DIaz));
    BWaz=zeros(size(DIaz));
    BW3az=zeros(size(DIaz));
    BW15az=zeros(size(DIaz));
    Bandwidthaz = nan(2,1);
    if strcmp(Config.ArrayMan.SteerSpace,'3D')  
        Resp3d=zeros(nfreqs,nlook);
        DI3d=zeros(nfreqs,1);
        WNG=zeros(size(DIaz));
        Contrast3d=zeros(size(DIaz));
        ContrastRange3d=zeros(size(DIaz));
    end
end        
if strcmp(Config.ArrayMan.SteerSpace,'2Del')
    Respdirel=zeros(nfreqs,nel);
    Respshiftel=zeros(nfreqs,nel);
    Respel=zeros(nfreqs,nel);
    RespOffTarEl = nan(nfreqs,1);

    DIel=zeros(nfreqs,1);
    WNG=zeros(size(DIel));
    Contrastel=zeros(size(DIel));
    ContrastRangeel=zeros(size(DIel));
    SLel=zeros(size(DIel));
    BWel=zeros(size(DIel));
    BW3el=zeros(size(DIel));
    BW15el=zeros(size(DIel));
    Bandwidthel = nan(2,1);
end

w=BFWeights.fd;

GLTh = 0.1;  %Grating Lobe Threshold from the main lobe in dB
SSLTh = 0;    %Sidelobe Suppression Threshold in dB between SSL along angle vs that along f
deltafindTh = 0;
kmaxSSL = round(nfreqs/3);    %Maximum frequency index to check that SSL is also a local maximum along f (only for low frequencies)
nmaxsidelobes = 5;
nSidelobes = 20;

if strcmp(Config.ArrayMan.SteerSpace,'3D')
    AngResAz = ArrayMan(1).LookAng.Az(2)-ArrayMan(1).LookAng.Az(1);
    AngResEl = ArrayMan(1).LookAng.El(2)-ArrayMan(1).LookAng.El(1);
elseif strcmp(Config.ArrayMan.SteerSpace,'2Daz')
    AngResAz = ArrayMan(1).LookAng.Az(2)-ArrayMan(1).LookAng.Az(1);
    AngResEl = 1;
elseif strcmp(Config.ArrayMan.SteerSpace,'2Del')
    AngResEl = ArrayMan(1).LookAng.El(2)-ArrayMan(1).LookAng.El(1);
    AngResAz = 1;
end

%running loop backwards for allocating memory in first iteration
for inoise=size(ArrayMan,1):-1:1
    for itest=size(ArrayMan,2):-1:1
    if inoise==1; itest=1; end
    
    a=ArrayMan(inoise,itest).a;
    Rvv=ArrayMan(inoise,itest).Rvv;
    
    GLfcount = 0;
    GratingLobesAz = [];
    GratingLobesEl = [];
    GratingLobesAzf = [];
    GratingLobesElf = [];
    sidelobesAz = nan(nfreqs,nSidelobes*2);
    sidelobeAzind = zeros(nfreqs,1);
    nAzsidelobes = zeros(nfreqs,1);
    sidelobeElind = zeros(nfreqs,1);
    nElsidelobes = zeros(nfreqs,1);
    
    for k=1:nfreqs
%         [~,targetAzind]=min(abs(deg2rad(targetAz)-ArrayMan(1).LookAng.Az));
%         [~,targetElind]=min(abs(deg2rad(targetEl)-ArrayMan(1).LookAng.El));
%         
        atmp=a(:,:,k); 
        
        %%% for azimuth only and 3D space
        if strcmp(Config.ArrayMan.SteerSpace,'2Daz')
        aaz=atmp(:,iaz); %#ok<*FNDSB>
        aaz_target=atmp(:,iazDirTarget);
        if doInterf; aaz_int=atmp(:,iazInt); end

        Respdiraz(k,:)=abs(w(:,k)'*aaz);
        Respshiftaz(k,:) = cat(2,Respdiraz(k,iaz0:end),Respdiraz(k,1:iaz0-1));
        [~,peakdirinds]=findpeaks(abs(Respdiraz(k,:)),'sortstr','descend','npeaks',20);
        [~,peakshiftinds]=findpeaks(abs(Respshiftaz(k,:)),'sortstr','descend','npeaks',20);
        [~,troughdirinds]=findpeaks(-db(Respdiraz(k,:)));
        [~,troughshiftinds]=findpeaks(-db(Respshiftaz(k,:)));
        if iazDirTarget < iaz270 || iazDirTarget > iaz90
            Respaz(k,:) = Respshiftaz(k,:);
            %peakinds = peakwrapinds+Az0ind-1;
%             shift=0;
%             for ipeak=1:length(peakinds)
%                 if peakinds(ipeak)>Az180ind
%                     shift=1;
%                 end
%             end
%             if shift
%                 peakinds = peakwrapinds-Az0ind-1;
%             end
            shift=1;    %used later to take the shifted response
            if iazDirTarget < iaz270
                iazTarget = iazDirTarget +iaz0+1;
            else
                iazTarget = iazDirTarget -iaz0+1;
            end
        else
            iazTarget = iazDirTarget;
            Respaz(k,:) = Respdiraz(k,:);
            % this attempts to ignore aliasing effects from grating lobes
            shift=0;
        end
        peakinds = zeros(1,length(peakshiftinds));
        if ~isempty(peakshiftinds)
            for ip = 1:length(peakshiftinds)
                if peakshiftinds(ip) <= iaz0+1  %if it is between 0 and 180
                    peakinds(ip) = peakshiftinds(ip) +iaz0-1;
                else    %if it is between 180 and 360
                    peakinds(ip) = peakshiftinds(ip) -iaz0-1;
                end
            end
        end
        troughinds = zeros(1,length(troughshiftinds));
        if ~isempty(troughshiftinds)
            for ip = 1:length(troughshiftinds)
                if troughshiftinds(ip) <= iaz0+1  %if it is between 0 and 180
                    troughinds(ip) = troughshiftinds(ip) +iaz0-1;
                else    %if it is between 180 and 360
                    troughinds(ip) = troughshiftinds(ip) -iaz0-1;
                end
            end
        end
        
        peakinds = unique(cat(2,peakinds, peakdirinds));
        peakinds = peakinds(peakinds~=0);
        [~,minpeakdiffind] = min(abs(iDirTarget-peakinds));
        mainPeakind = peakinds(minpeakdiffind);
        mainPeak = db(Respaz(k,mainPeakind));
        
        troughinds = unique(cat(2,troughinds, troughdirinds));
        troughinds = troughinds(troughinds~=0);
        
        if ~isempty(mainPeak)
            
        %calculating potential grating lobes (assumed to be a maximum of 2)
        GratingLobesind = find(db(Respdiraz(k,peakinds~=mainPeakind))>mainPeak-GLTh,2);

        if ~isempty(GratingLobesind)
            if ~exist('GratingLobesAz','var')
                GratingLobesAz = zeros(nfreqs,2);
                GratingLobesAzf = zeros(nfreqs,1);
            end
            GLfcount = GLfcount + 1;
            if shift
                GratingLobesind = GratingLobesind +iaz0-1;
            end
            GratingLobesAz(GLfcount,1:length(GratingLobesind)) = (GratingLobesind-1)*AngResAz + ArrayMan(1).LookAng.Az(1);
            GratingLobesAzf(GLfcount) = ArrayMan(1).f(k);
        end
        
        
        % 3dB beamwidth
            %BeamwidthEstimateAzind=abs(find(Respaz(k,1:targetAzind)<=mainPeak-3,1,'last')-(targetAzind+find(Respaz(k,targetAzind:end)<=mainPeak-3,1,'first') - 1));
            Beamwidth3EstimateAzind=abs(find(db(Respaz(k,1:iazTarget))<=db(Respaz(k,iazTarget))-3,1,'last')-(iazTarget+find(db(Respaz(k,iazTarget:end))<=db(Respaz(k,iazTarget))-3,1,'first') - 1));
            Beamwidth3EstimateAz = Beamwidth3EstimateAzind*AngResAz;
            Beamwidth15EstimateAzind=abs(find(db(Respaz(k,1:iazTarget))<=db(Respaz(k,iazTarget))-15,1,'last')-(iazTarget+find(db(Respaz(k,iazTarget:end))<=db(Respaz(k,iazTarget))-15,1,'first') - 1));
            Beamwidth15EstimateAz = Beamwidth15EstimateAzind*AngResAz;
            null1ind = min(troughinds(iazTarget<troughinds));
            null2ind = max(troughinds(iazTarget>troughinds));
            BeamwidthEstimateAz = (null1ind-null2ind)*AngResAz;
        %end
        
        if isempty(Beamwidth3EstimateAz)
            BW3az(k)=360;
        else
            BW3az(k)=Beamwidth3EstimateAz;
        end
        if isempty(Beamwidth15EstimateAz)
            BW15az(k)=360;
        else
            BW15az(k)=Beamwidth15EstimateAz;
        end
        if isempty(BeamwidthEstimateAz)
            BWaz(k)=360;
        else
            BWaz(k)=BeamwidthEstimateAz;
        end
        
        % Side-lobe suppression (main lobe level - highest side lobe level)
        if isempty(peakinds)||numel(peakinds)<2 %modified to 2 for asymmetric beampatterns like MVDR
            SLaz(k)=0;
        else
            %the position of the sidelobe is after all the grating lobes +
            %the main peak + 1
            %sidelobeinds = peakinds(length(GratingLobesind) + 2 :end);            
%             if isempty(sidelobeinds)||numel(sidelobeinds)<3
%                 SLaz(k)=0;
%             else
%                 %SLaz(k)=abs(Respaz(k,mainPeakind)-max(Respaz(k,sidelobeinds(1))));
%                 SLaz(k)=abs(Respaz(k,targetAzind)-max(Respaz(k,sidelobeinds(1))));
%             end
            
            %testing with grating lobes included as sidelobes
            sidelobeinds = peakinds(peakinds~=mainPeakind);
            [sidelobesAz(k,1:length(sidelobeinds)),sidelobeorder] = sort(abs(Respdiraz(k,sidelobeinds)),'descend');
            sidelobeinds = sidelobeinds(sidelobeorder);
            nAzsidelobes(k) = length(sidelobeinds);
            sidelobeAzind(k) = sidelobeinds(1);
            SLaz(k)=abs(db(Respaz(k,iTarget))-db(Respdiraz(k,sidelobeinds(1))));
            SLazInd(k) = sidelobeinds(1);
            RespOffTarAz(k) = sum(abs(Respaz(k,1:null2ind)).^2) + sum(10.^(Respaz(k,null1ind:end)/10).^2);
            
        end
        
        else
            BWaz(k)=360;
            BW3az(k)=360;
            BW15az(k)=360;
            SLaz(k)=0;
        end
        
        % Array manifold vector steered at the target/interefer(s) direction
        % DIaz only when 2D azimuth, otherwise it will be DI3d (below)
        Rvvtmp=squeeze(Rvv(:,:,k));
        DIaz(k)=10*log10(abs(w(:,k)'*aaz_target)^2./real(w(:,k)'*Rvvtmp*w(:,k)));
        WNG(k)=getcurrWNG(w(:,k),aaz_target);
        %WNG(k)=10*log10(abs(w(:,k)'*a_target)^2./(w(:,k)'*w(:,k)));
        
        %response at each of the interferers
        Respintaz=zeros(1,nints);
        for kint=1:nints
            Respintaz(kint)=abs(w(:,k)'*squeeze(aaz_int(:,kint)))^2;
        end
        Respintmeanaz=mean(Respintaz);
        %Mean Acoustic contrast between target and interferers
        Contrastaz(k)=10*log10(abs(w(:,k)'*aaz_target)^2/Respintmeanaz);
        ContrastRangeaz(k) = 10*log10(sum(Respaz(k,iazRangemin:iazRangemax).^2,2)./...
            sum(cat(2,Respaz(k,1:iazRangemin-1).^2,Respaz(k,iazRangemax+1:end).^2),2));
        if shift
            Respaz(k,:)=Respdiraz(k,:);
        end
        
        %%% Response and metrics for elevation and 3D
        elseif strcmp(Config.ArrayMan.SteerSpace,'2Del')
        
        ael=squeeze(atmp(:,iel));
        ael_target=atmp(:,ielDirTarget);
        if doInterf; ael_int=atmp(:,ielInt); end
        Respdirel(k,:)=w(:,k)'*ael;
        Respshiftel(k,:) = cat(2,Respdirel(k,iel0:end),Respdirel(k,1:iel0-1));
        [~,peakdirinds]=findpeaks(abs(Respdirel(k,:)),'sortstr','descend','npeaks',20);
        [~,peakshiftinds]=findpeaks(abs(Respshiftel(k,:)),'sortstr','descend','npeaks',20);
        [~,troughdirinds]=findpeaks(-db(Respdirel(k,:)));
        [~,troughshiftinds]=findpeaks(-db(Respshiftel(k,:)));
        if ielDirTarget < ielminus90 || ielDirTarget > iel90
            Respel(k,:) = Respshiftel(k,:);
            %peakinds = peakwrapinds+El0ind-1;
%             shift=0;
%             for ipeak=1:length(peakinds)
%                 if peakinds(ipeak)>El180ind
%                     shift=1;
%                 end
%             end
%             if shift
%                 peakinds = peakwrapinds-El0ind-1;
%             end
            shift=1;    %used later to take the shifted response
            if ielDirTarget < ielminus90
                ielTarget = ielDirTarget +iel0+1;
            else
                ielTarget = ielDirTarget -iel0+1;
            end
        else
            ielTarget = ielDirTarget;
            Respel(k,:) = Respdirel(k,:);
            % this attempts to ignore aliasing effects from grating lobes
            shift=0;
        end
        peakinds = zeros(1,length(peakshiftinds));
        if ~isempty(peakshiftinds)
            for ip = 1:length(peakshiftinds)
                if peakshiftinds(ip) <= iel0+1  %if it is between 0 and 180
                    peakinds(ip) = peakshiftinds(ip) +iel0-1;
                else    %if it is between 180 and 360
                    peakinds(ip) = peakshiftinds(ip) -iel0-1;
                end
            end
        end
        troughinds = zeros(1,length(troughshiftinds));
        if ~isempty(troughshiftinds)
            for ip = 1:length(troughshiftinds)
                if troughshiftinds(ip) <= iel0+1  %if it is between 0 and 180
                    troughinds(ip) = troughshiftinds(ip) +iel0-1;
                else    %if it is between 180 and 360
                    troughinds(ip) = troughshiftinds(ip) -iel0-1;
                end
            end
        end    
            
        if max(el)<=90 && min(el)>=-90
        if ielDirTarget < ielminus90
            ielTarget = ielDirTarget +iel0+1;
            %taking the bottom half, mirroring and concatenating it with
            %half we already have
            %since we are missing elevation -90 we assume this to be the
            %same as -90 + AngResEl for now
            Respel(k,:) = cat(2,Respel(k,iel0-1:-1:1),Respel(k,1),Respel(k,1:iel0));
            [~,peakinds]=findpeaks(Respel(k,:),'sortstr','descend','npeaks',20);
            shift=1;    %used later to take the shifted response
        elseif ielDirTarget > iel90
            ielTarget = ielDirTarget -iel0+1;
            Respel(k,:) = cat(2,Respel(k,iel0:end),Respel(k,end-1:-1:iel0+1));
            [~,peakinds]=findpeaks(Respel(k,:),'sortstr','descend','npeaks',20);
            shift=1;    %used later to take the shifted response
        else
            % this attempts to ignore aliasing effects from grating lobes
            [~,peakinds]=findpeaks(Respel(k,:),'sortstr','descend','npeaks',20);
            shift=0;
        end
        end
        
        peakinds = unique(cat(2,peakinds, peakdirinds));
        peakinds = peakinds(peakinds~=0);
        [~,minpeakdiffind] = min(abs(ielDirTarget-peakinds));
        mainPeakind = peakinds(minpeakdiffind);
        mainPeak = db(Respel(k,mainPeakind));
        
        troughinds = unique(cat(2,troughinds, troughdirinds));
        troughinds = troughinds(troughinds~=0);
        
        if ~isempty(mainPeak)
        %calculating potential grating lobes (assumed to be a maximum of 2)
        GratingLobesind = find(db(Respdirel(k,peakinds~=mainPeakind))>mainPeak-GLTh,2);
        
        if ~isempty(GratingLobesind)
            if ~exist('GratingLobesEl','var')
                GratingLobesEl = zeros(nfreqs,2);
                GratingLobesElf = zeros(nfreqs,1);
            end
            GLfcount = GLfcount + 1;
            if shift
                GratingLobesind = GratingLobesind +iel0-1;
            end
            GratingLobesEl(GLfcount,1:length(GratingLobesind)) = (GratingLobesind-1)*AngResEl + ArrayMan(1).LookAng.El(1);
            GratingLobesElf(GLfcount) = ArrayMan(1).f(k);
        end
        
        % 3dB beamwidth
            %BeamwidthEstimateElind=abs(find(Respel(k,1:targetElind)<=mainPeak-3,1,'last')-(targetElind+find(Respel(k,targetElind:end)<=mainPeak-3,1,'first') - 1));
            Beamwidth3EstimateElind=abs(find(db(Respel(k,1:ielTarget))<=db(Respel(k,ielTarget))-3,1,'last')-(ielTarget+find(db(Respel(k,ielTarget:end))<=db(Respel(k,ielTarget))-3,1,'first') - 1));
            Beamwidth3EstimateEl = Beamwidth3EstimateElind*AngResEl;
            Beamwidth15EstimateElind=abs(find(db(Respel(k,1:ielTarget))<=db(Respel(k,ielTarget))-15,1,'last')-(ielTarget+find(db(Respel(k,ielTarget:end))<=db(Respel(k,ielTarget))-15,1,'first') - 1));
            Beamwidth15EstimateEl = Beamwidth15EstimateElind*AngResEl;
            null1ind = min(troughinds(ielTarget<troughinds));
            null2ind = max(troughinds(ielTarget>troughinds));
            BeamwidthEstimateEl = (null1ind-null2ind)*AngResEl;
        %end
        
        if isempty(Beamwidth3EstimateEl)
            BW3el(k)=360;
        else
            BW3el(k)=Beamwidth3EstimateEl;
        end
        if isempty(Beamwidth15EstimateEl)
            BW15el(k)=360;
        else
            BW15el(k)=Beamwidth15EstimateEl;
        end
        if isempty(BeamwidthEstimateEl)
            BWel(k)=360;
        else
            BWel(k)=BeamwidthEstimateEl;
        end
        
        % Side-lobe suppression (main lobe level - highest side lobe level)
        if isempty(peakinds)||numel(peakinds)<2 %modified to 2 for asymmetric beampatterns like MVDR
            SLel(k)=0;
        else
            %the position of the sidelobe is after all the grating lobes +
            %the main peak + 1
            %sidelobeinds = peakinds(length(GratingLobesind) + 2 :end);            
%             if isempty(sidelobeinds)||numel(sidelobeinds)<3
%                 SLel(k)=0;
%             else
%                 %SLel(k)=abs(Respel(k,mainPeakind)-max(Respel(k,sidelobeinds(1))));
%                 SLel(k)=abs(Respel(k,targetElind)-max(Respel(k,sidelobeinds(1))));
%             end
            
            %testing with grating lobes included as sidelobes            
            sidelobeinds = peakinds(peakinds~=mainPeakind);
            [sidelobesEl(k,1:length(sidelobeinds)),sidelobeorder] = sort(db(Respdirel(k,sidelobeinds)),'descend');
            sidelobeinds = sidelobeinds(sidelobeorder);
            nElsidelobes(k) = length(sidelobeinds);
            sidelobeElind(k) = sidelobeinds(1);
            SLel(k)=abs(db(Respel(k,iTarget))-db(Respdirel(k,sidelobeinds(1))));
            SLelInd(k) = sidelobeinds(1);
            RespOffTarEl(k) = sum(abs(Respel(k,1:null2ind)).^2) + sum(abs(Respel(k,null1ind:end)).^2);
        end
        
        else
            BWel(k)=360;
            BW3el(k)=360;
            BW15el(k)=360;
            SLel(k)=0;
        end
        

        % Array manifold vector steered at the target/interefer(s) direction
        % DIel only when 2D elevation, otherwise it will be DI3d (below)
        if strcmp(Config.ArrayMan.SteerSpace,'2Del')      
            Rvvtmp=squeeze(Rvv(:,:,k));
            DIel(k)=10*log10(abs(w(:,k)'*ael_target)^2./real(w(:,k)'*Rvvtmp*w(:,k)));
        end
        WNG(k)=getcurrWNG(w(:,k),ael_target);
        
        %response at each of the interferers
        Respintel=zeros(1,nints);
        for kint=1:nints
            Respintel(kint)=abs(w(:,k)'*ael_int(:,kint))^2;
        end
        Respintmeanel=mean(Respintel);
        %Mean Acoustic contrast between target and interferers
        Contrastel(k)=10*log10(abs(w(:,k)'*ael_target)^2/Respintmeanel);
        ContrastRangeel(k) = 10*log10(sum(Respel(k,ielRangemin:ielRangemax).^2,2,'omitnan')./...
            sum(cat(2,Respel(k,1:ielRangemin-1).^2,Respel(k,ielRangemax+1:end).^2),2,'omitnan'));
        
        if shift
            Respel(k,:)=Respdirel(k,:);
        end
        
        %%% for 3D space
        elseif strcmp(Config.ArrayMan.SteerSpace,'3D')   
        a_target=atmp(:,iTarget);
        a_int=atmp(:,iInt);
        Rvvtmp=squeeze(Rvv(:,:,k));
        Resp3d(k,:)=w(:,k)'*atmp;
        DI3d(k)=10*log10(abs(w(:,k)'*a_target)^2./real(w(:,k)'*Rvvtmp*w(:,k)));
        WNG(k)=getcurrWNG(w(:,k),a_target);

        %response at each of the interferers
        Respint=zeros(1,nints);
        for kint=1:nints
            Respint(kint)=abs(w(:,k)'*squeeze(a_int(:,kint)))^2;
        end
        Respintmean=mean(Respint);
        %Mean Acoustic contrast between target and interferers
        Contrast3d(k)=10*log10(abs(w(:,k)'*a_target)^2/Respintmean);
        ContrastRange3d(k) = 20*log10(sum(abs(Resp3d(k,iRangeOn)),'omitnan')./sum(abs(Resp3d(k,iRangeOff)),'omitnan'));
%         ContrastRange3D(k) = 20*log10(sum(sum(10.^(Resp3D(k,elRangeminind:elRangemaxind,azRangeminind:azRangemaxind)/20),2),3)./...
%             sum(sum(10.^(cat(2,Resp3D(k,1:elRangeminind-1,1:azRangeminind-1)/20,Respaz(k,elRangemaxind+1:end,azRangemaxind+1:end)/20)),2),3));
        end
    end
    
    %need to run another loop as Resp was not calculated for all freq in
    %previous loop to check that SSL is also a local maximum along freq
    if strcmp(Config.ArrayMan.SteerSpace,'2Daz')      
    for k=1:kmaxSSL
        %if more than nmaxsidelobes then no need to check local maximum
        %along f as uncertainty is only with one sidelobe on either side of
        %mainlobe
%         if nsidelobes(k)>=nmaxsidelobes
%             break
%         end
        if SLaz(k)~=0
            %checking that the sidelobe is a local maximum for all frequencies,
            %otherwise should not be considered
            [~,peakfinds]=findpeaks(abs(Respdiraz(1:kmaxSSL,sidelobeAzind(k))),'npeaks',5);
            [deltafind,peakfind] = min(abs(peakfinds-k));
            
            if ~isempty(peakfind) && length(peakfinds)>1
                %threshold set as a quarter of the distance between two peaks
                %deltafindTh = (peakfinds(2) - peakfinds(1))/4;
                if db(Respdiraz(k,sidelobeAzind(k)))-db(Respdiraz(peakfinds(peakfind),sidelobeAzind(k)))>SSLTh || deltafind>deltafindTh
                    SLaz(k)=0;
                else 
                    break
                end
            end
        end
    end
    end
    
    %bandwidth
    if strcmp(Config.ArrayMan.SteerSpace,'2Daz')
        for iphi=1:naz
            [~,ipeaks] = findpeaks(-db(Respaz(:,iphi)),'NPeaks',30);
            ipeaks = sort(ipeaks);
            if isempty(ipeaks); ipeaks = length(ArrayMan(1).f); end
            freqminind(iphi) = ipeaks(1);
        end
        freqminind = min(freqminind((1:naz)~=iazTarget));
        BWaz(1:freqminind-1) = 360;
        if strcmp(Config.Array.Geo,'linear_bs')
            LFind = find(BWaz<180,1,'first');
        else
            LFind = find(BWaz<360,1,'first');
        end
%         if ~isempty(GratingLobesAzf)
%             GratingLobesAzf = GratingLobesAzf(GratingLobesAzf> ArrayMan(1).f(LFind));
%             if ~isempty(GratingLobesAzf)
%                 [~,HFind] = min(abs(min(GratingLobesAzf) - ArrayMan(1).f));
%             end
%         end
%         if isempty(GratingLobesAzf)
%             ifreq = find(sum(isnan(sidelobesAz),2)~=size(sidelobesAz,2),1);
%             sumSideLobes = sum(10.^(sidelobesAz/10),2,'omitnan');
%             navg = 20;
%             
%             for n=1:length(sumSideLobes)-navg
%                 sumSideLobeSmooth(n) = median(sumSideLobes(n:navg+n-1));
%             end
%             [~,HFind] = findpeaks(sumSideLobeSmooth,'MinPeakHeight',1.2.*sumSideLobeSmooth(1),'NPeaks',20,'SortStr','none');
%             HFind = HFind(1) + round(navg/2) + ifreq-1;
%             %HFind = find(diff()>0.1*sum(sidelobes(2:end),2));
%             %[~,HFind] = max(diff(abs(SLazInd-targetDirAzind)));
%         end
%         if HFind==1; HFind=length(ArrayMan(1).f); end
        RespOffTarAz = RespOffTarAz ./ abs(Respaz(:,iazTarget)).^2;
        RespOffTarAz = RespOffTarAz/max(RespOffTarAz);
        [OffTarPeaks,OffTarPeakPos] = findpeaks(RespOffTarAz,'MinPeakProminence',0.01);
        [~,OffTarPromPeakPos] = findpeaks(RespOffTarAz,'MinPeakProminence',0.05);
        HFind=length(ArrayMan(1).f);
        for ipeak = 1:length(OffTarPromPeakPos)
            if sum(OffTarPromPeakPos(ipeak)==OffTarPeakPos([false;OffTarPeaks(2:end)>OffTarPeaks(1:end-1)]))>0
                HFind = OffTarPromPeakPos(ipeak);
                break;
            end
        end
        %HFind = LFind -1 + find(10*log10(abs(RespOffTarAz(LFind:end)))>-10,1,'first');
        %if isempty(HFind); HFind=length(ArrayMan(1).f); end
        if ~isempty(LFind); Bandwidthaz(1) = ArrayMan(1).f(LFind); end
        if ~isempty(HFind); Bandwidthaz(2) = ArrayMan(1).f(HFind); end

    elseif strcmp(Config.ArrayMan.SteerSpace,'2Del') 
        for itheta=1:nel
            [~,ipeaks] = findpeaks(-db(Respel(:,itheta)),'NPeaks',30);
            ipeaks = sort(ipeaks);
            if isempty(ipeaks); ipeaks = length(ArrayMan(1).f); end
            freqminind(itheta) = ipeaks(1);
        end
        freqminind = min(freqminind((1:nel)~=iTarget));
        BWel(1:freqminind-1) = 360;
        if strcmp(Config.Array.Geo,'linear_bs')
            LFind = find(BWel<180,1,'first');
        else
            LFind = find(BWel<360,1,'first');
        end
        RespOffTarEl = RespOffTarEl ./ abs(Respel(:,iTarget)).^2;
        RespOffTarEl = RespOffTarEl/max(RespOffTarEl);
        [OffTarPeaks,OffTarPeakPos] = findpeaks(RespOffTarEl,'MinPeakProminence',0.01);
        [~,OffTarPromPeakPos] = findpeaks(RespOffTarEl,'MinPeakProminence',0.05);
        HFind=length(ArrayMan(1).f);
        for ipeak = 1:length(OffTarPromPeakPos)
            if sum(OffTarPromPeakPos(ipeak)==OffTarPeakPos([false;OffTarPeaks(2:end)>OffTarPeaks(1:end-1)]))>0
                HFind = OffTarPromPeakPos(ipeak);
                break;
            end
        end
        %HFind = LFind -1 + find(10*log10(abs(RespOffTarAz(LFind:end)))>-10,1,'first');
        %if isempty(HFind); HFind=length(ArrayMan(1).f); end
        if ~isempty(LFind); Bandwidthel(1) = ArrayMan(1).f(LFind); end
        if ~isempty(HFind); Bandwidthel(2) = ArrayMan(1).f(HFind); end
    end
    
    %passing variables to structure
    if strcmp(Config.ArrayMan.SteerSpace,'2Daz')
        Eval(inoise,itest).DIaz=DIaz;
        Eval(inoise,itest).BWaz=BWaz;
        Eval(inoise,itest).BW3az=BW3az;
        Eval(inoise,itest).BW15az=BW15az;
        Eval(inoise,itest).SLaz=SLaz;
        Eval(inoise,itest).Contrastaz=Contrastaz;
        Eval(inoise,itest).ContrastRangeaz = ContrastRangeaz;
        Eval(inoise,itest).WNG=WNG;
        Eval(inoise,itest).Respaz=Respaz;
        Eval(inoise,itest).GL.Az=GratingLobesAz(GratingLobesAzf~=0);
        Eval(inoise,itest).GL.Azf=GratingLobesAzf(GratingLobesAzf~=0);
        Eval(inoise,itest).FRangeaz = Bandwidthaz;
        Eval(inoise,itest).FRespaz = Respaz(:,iTarget);
 
    elseif strcmp(Config.ArrayMan.SteerSpace,'2Del') 
        Eval(inoise,itest).DIel=DIel;
        Eval(inoise,itest).BWel=BWel;
        Eval(inoise,itest).BW3el=BW3el;
        Eval(inoise,itest).BW15el=BW15el;
        Eval(inoise,itest).SLel=SLel;
        Eval(inoise,itest).Contrastel=Contrastel;
        Eval(inoise,itest).ContrastRangeel = ContrastRangeel;
        Eval(inoise,itest).Respel=Respel;
        Eval(inoise,itest).GL.El=GratingLobesEl(GratingLobesElf~=0);
        Eval(inoise,itest).GL.Elf=GratingLobesElf(GratingLobesElf~=0);
        Eval(inoise,itest).FRangeel = Bandwidthel;
        Eval(inoise,itest).FRespel = Respel(:,iTarget);
        Eval(inoise,itest).WNG=WNG;
    elseif strcmp(Config.ArrayMan.SteerSpace,'3D') 
        Eval(inoise,itest).DI3d=DI3d;
        Eval(inoise,itest).WNG=WNG;
        Eval(inoise,itest).Resp3d=Resp3d;
        Eval(inoise,itest).Contrast3d=Contrast3d;
        Eval(inoise,itest).ContrastRange3d = ContrastRange3d;
    end
    
    end
    
end

end %function