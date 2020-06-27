function DirectionalResponseMeasuresFreqAvg(ConfigAllArrays,BFEvalAllArrays, FigConfig)
% function DirectionalResponseMeasuresFreqAvg(ConfigAllArrays,BFEvalAllArrays, FigConfig)
% Function that plots performance metrics from the beamforming analysis as
% the average value across frequency.
%
% input arguments:
%   ConfigAllArrays: cell array containing all the configuration settings 
%       for each microphone array, beamforming method, etc.
%   BFEvalAllArrays: struct containing containing several performance
%       metrics for each array, beamforming method, etc. Note 
%       BFEvalAllArrays is a struct from the output of 
%       rearrange_directional_response_measures.m. 
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.


f=ConfigAllArrays{1,1,1,1}.Filter.FreqVector;
freqs = length(f);
foctaves = [32 63 125 250 500 1000 2000 4000 8000 16000];
fbandcutoff = [foctaves.*2^(-1/2) foctaves(end).*2^(1/2)];
%f=linspace(fbandcutoff(1),fbandcutoff(end),1000);  %test to check that
%with linear frequency vector starting and ending at fbandcutoff the number
%of frequencies per octave band is constant.
ifirstoctave = find(foctaves>=f(1),1,'first');
ilastoctave = find(foctaves<=f(end),1,'last');
foctaves = foctaves(ifirstoctave:ilastoctave);
iband = assign_freq_to_bands(f,fbandcutoff);

flog = logspace(log10(f(1)),log10(f(end)),20*length(foctaves));
flog_prime = zeros(1,length(flog));
iflog_prime = zeros(1,length(flog));
for iflog=1:length(flog)
    [flog_prime(iflog),iflog_prime(iflog)] = min(abs(f-flog(iflog)));
end
ibandlog = assign_freq_to_bands(flog,fbandcutoff);
countbandlog = zeros(length(foctaves),1);
for iflog=1:length(flog)
    countbandlog(ibandlog(iflog)) = countbandlog(ibandlog(iflog)) +1;
end
    
narray_geo = size(ConfigAllArrays,1);
nM = size(ConfigAllArrays,2);
nr = size(ConfigAllArrays,3);
nepsilons = size(ConfigAllArrays,4);
nmethods = size(ConfigAllArrays{1,1,1}.Filter.Method,2);
nsources = size(ConfigAllArrays{1,1,1}.Scene.Targets,1);
nnoises = size(BFEvalAllArrays,6);
%nepsilons = size(BFEvalAllArrays,7);

source1 = 1;
sourcen = 1;
epsilon1 = 1;
epsilonn = 1;
noise1 = 1;
noisen = nnoises;

array_geo = cell(1,narray_geo);
array_geo_disp = cell(1,narray_geo);

M = zeros(1,nM);
r = zeros(1,nr);
%radius estimated from half the aperture of linear array
for ir=1:nr
    r(ir)= round(ConfigAllArrays{1,1,ir}.Array.Aperture/2,2);
end
rmm = r*1000;
fmin = zeros(size(ConfigAllArrays));
fmax = zeros(size(ConfigAllArrays));
ifmin = zeros(size(ConfigAllArrays));
ifmax = zeros(size(ConfigAllArrays));
for iarray_geo=1:narray_geo
    array_geo{iarray_geo} = ConfigAllArrays{iarray_geo,1,1}.Array.Geo;
    array_geo_disp{iarray_geo} = ConfigAllArrays{iarray_geo,1,1}.Array.GeoDisp;
    for iM=1:nM
        for ir=1:nr
            fmin(iarray_geo,iM,ir) = ConfigAllArrays{iarray_geo,iM,ir}.ArrayMan.c/ConfigAllArrays{iarray_geo,iM,ir}.Array.Aperture;
            fmax(iarray_geo,iM,ir) = ConfigAllArrays{iarray_geo,iM,ir}.Array.falias;
            [~,ifmin(iarray_geo,iM,ir)] = min(abs(fmin(iarray_geo,iM,ir)-f));
            [~,ifmax(iarray_geo,iM,ir)] = min(abs(fmax(iarray_geo,iM,ir)-f));
            ConfigAllArrays{iarray_geo,iM,ir}.fcuton = fmin(iarray_geo,iM,ir);
        end
    end
end

for iM=1:nM
    M(iM) = size(ConfigAllArrays{1,iM,1}.Array.MicPos,1);
end


methodStr = ConfigAllArrays{1,1,1}.Filter.Method;
doaz=0;
doel=0;
do3d=0;
if isfield(BFEvalAllArrays,'BWaz')
    DIazbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    BWazbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    BW15azbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    SLazbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    WNGbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    ACazbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    
%     if isfield(BFEvalAllArrays,'FRespaz')
%         doFR=1;
%         FRangeaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,2);
%         FRespaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
%     end
    if isfield(BFEvalAllArrays,'ContrastRangeaz')
        ACRangeazbands = nan(length(foctaves),freqs);
        doACRange=1;
    end
    doaz=1;
end

if isfield(BFEvalAllArrays,'BWel')

    DIelbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    BWelbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    BW15elbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    SLelbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    WNGbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    ACelbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
    
%     if isfield(BFEvalAllArrays,'FRespaz')
%         doFR=1;
%         FRangeaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,2);
%         FRespaz = zeros(narray_geo,nM,nr,nmethods,nsources,nnoises,nepsilons,freqs);
%     end
    if isfield(BFEvalAllArrays,'ContrastRangeaz')
        ACRangeelbands = nan(length(foctaves),freqs);
        doACRange=1;
    end
    doel=1;
end
if isfield(BFEvalAllArrays,'DI3d')
    DI3d = nan(length(foctaves),freqs);
    do3d=1;
end    

DIelbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
DI3dbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
BWazbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
BWelbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
BW15azbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
BW15elbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
SLazbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
SLelbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
WNGbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
ACazbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)
ACelbands = nan(length(foctaves),freqs);    %2D matrix (ibands x however many freq/band)

for iarray_geo=1:narray_geo
    for iM=1:nM
        for ir=1:nr
            for imethod=1:nmethods
                for isource=1:nsources
                    for inoise=1:nnoises
                        for ie=1:nepsilons
                        count = zeros(length(foctaves),1);    %position of each band for each freq added
                        for ifreq=ifmin(iarray_geo,iM,ir):ifmax(iarray_geo,iM,ir)
                            count(iband(ifreq)) = count(iband(ifreq)) + 1;
                            if doaz
                                DIazbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.DIaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                                BWazbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.BWaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                                BW15azbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.BW15az(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                                SLazbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.SLaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                                WNGbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                                ACazbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.Contrastaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);                          
                                if doACRange; ACRangeazbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.ContrastRangeaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq); end                            
                            end
                            if doel
                                DIelbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.DIel(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                                BWelbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.BWel(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                                BW15elbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.BW15el(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                                SLelbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.SLel(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                                WNGbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                                ACelbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.Contrastaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);                          
                                if doACRange; ACRangeelbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.ContrastRangeel(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq); end                            
                            end
                            if do3d
                                 DI3dbands(iband(ifreq),count(iband(ifreq))) = BFEvalAllArrays.DI3d(iarray_geo,iM,ir,imethod,isource,inoise,ie,ifreq);
                            end
                        end
                        if doaz
                            DIazfbands = mean(DIazbands,2,'omitnan');
                            BWazfbands = mean(BWazbands,2,'omitnan');
                            BW15azfbands = mean(BW15azbands,2,'omitnan');
                            SLazfbands = mean(SLazbands,2,'omitnan');
                            WNGfbands = mean(WNGbands,2,'omitnan');
                            ACazfbands = mean(ACazbands,2,'omitnan');
                            if doACRange
                                ACRangeazfbands = mean(ACRangeazbands,2,'omitnan');
                                ACRangeazfavg = mean(ACRangeazfbands,'omitnan');
                                BFEvalAllArrays.favg.ContrastRangeaz(iarray_geo,iM,ir,imethod,isource,inoise,ie) = ACRangeazfavg;
                            end
                            
                            DIazfavg = mean(DIazfbands,'omitnan');
                            BWazfavg = mean(BWazfbands,'omitnan');
                            BW15azfavg = mean(BW15azfbands,'omitnan');
                            SLazfavg = mean(SLazfbands,'omitnan');
                            WNGfavg = mean(WNGfbands,'omitnan');
                            ACazfavg = mean(ACazfbands,'omitnan');
                            
                            BFEvalAllArrays.favg.DIaz(iarray_geo,iM,ir,imethod,isource,inoise,ie) = DIazfavg;
                            BFEvalAllArrays.favg.BWaz(iarray_geo,iM,ir,imethod,isource,inoise,ie) = BWazfavg;
                            BFEvalAllArrays.favg.BW15az(iarray_geo,iM,ir,imethod,isource,inoise,ie) = BW15azfavg;
                            BFEvalAllArrays.favg.SLaz(iarray_geo,iM,ir,imethod,isource,inoise,ie) = SLazfavg;
                            BFEvalAllArrays.favg.WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie) = WNGfavg;
                            BFEvalAllArrays.favg.Contrastaz(iarray_geo,iM,ir,imethod,isource,inoise,ie) = ACazfavg;
                        end
                        
                        if doel
                            DIelfbands = mean(DIelbands,2,'omitnan');
                            BWelfbands = mean(BWelbands,2,'omitnan');
                            BW15elfbands = mean(BW15elbands,2,'omitnan');
                            SLelfbands = mean(SLelbands,2,'omitnan');
                            ACelfbands = mean(ACelbands,2,'omitnan');
                            WNGfbands = mean(WNGbands,2,'omitnan');
                            if doACRange
                                ACRangeelfbands = mean(ACRangeelbands,2,'omitnan'); 
                                ACRangeelfavg = mean(ACRangeelfbands,'omitnan');   
                                BFEvalAllArrays.favg.ContrastRangeel(iarray_geo,iM,ir,imethod,isource,inoise,ie) = ACRangeelfavg;
                            end
                            DIelfavg = mean(DIelfbands,'omitnan');
                            BWelfavg = mean(BWelfbands,'omitnan');
                            BW15elfavg = mean(BW15elfbands,'omitnan');
                            SLelfavg = mean(SLelfbands,'omitnan');
                            ACelfavg = mean(ACelfbands,'omitnan'); 
                            WNGfavg = mean(WNGfbands,'omitnan');
                            
                            BFEvalAllArrays.favg.DIel(iarray_geo,iM,ir,imethod,isource,inoise,ie) = DIelfavg;
                            BFEvalAllArrays.favg.BWel(iarray_geo,iM,ir,imethod,isource,inoise,ie) = BWelfavg;                        
                            BFEvalAllArrays.favg.BW15el(iarray_geo,iM,ir,imethod,isource,inoise,ie) = BW15elfavg;                        
                            BFEvalAllArrays.favg.SLel(iarray_geo,iM,ir,imethod,isource,inoise,ie) = SLelfavg;
                            BFEvalAllArrays.favg.Contrastel(iarray_geo,iM,ir,imethod,isource,inoise,ie) = ACelfavg;
                            BFEvalAllArrays.favg.WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie) = WNGfavg;
                        end
                            
                        if do3d
                            DI3dfbands = mean(DI3dbands,2,'omitnan');
                            DI3dfavg = mean(DI3dfbands,'omitnan');
                            BFEvalAllArrays.favg.DI3d(iarray_geo,iM,ir,imethod,isource,inoise,ie) = DI3dfavg;       
                        end

                        end
                    end
                end
            end
        end
    end
end



FigConfig.Fontsize = 14;
FigConfig.Xtickrotangle=35;
if ~isfield(FigConfig,'Save')
    FigConfig.Save = 0;
end

nmetrics = 7; 
ntotalmethods = 6;
do.metric = zeros(nmetrics);
do.bfmethod = zeros(ntotalmethods);
do.array = zeros(narray_geo);

for imetric=1:length(FigConfig.metrics)
    if strcmp(FigConfig.metrics{imetric}, 'all')
        do.metric(:) = 1; 
    elseif strcmp(FigConfig.metrics{imetric}, 'DI')
        do.metric(1) = 1;
    elseif strcmp(FigConfig.metrics{imetric}, 'WNG')
        do.metric(2) = 1;
    elseif strcmp(FigConfig.metrics{imetric}, 'BW') || strcmp(FigConfig.metrics{imetric}, 'BW3')
        do.metric(3) = 1;
    elseif strcmp(FigConfig.metrics{imetric}, 'BW15')
        do.metric(4) = 1;
    elseif strcmp(FigConfig.metrics{imetric}, 'SSL')    
        do.metric(5) = 1;
    elseif strcmp(FigConfig.metrics{imetric}, 'AC')    
        do.metric(6) = 1;
    elseif strcmp(FigConfig.metrics{imetric}, 'ACRange')
        do.metric(7) = 1;
    elseif strcmp(FigConfig.metrics{imetric}, 'FR')    
        do.metric(8) = 1;
    end
end 

for imethod=1:length(FigConfig.methods)
    if strcmp(FigConfig.methods{imethod}, 'all')
        do.bfmethod(:) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'ds')
        do.bfmethod(1) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'sda')
        do.bfmethod(2) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'mvdr')
        do.bfmethod(3) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'lcmv')    
        do.bfmethod(4) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'ls')    
        do.bfmethod(5) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'acc')    
        do.bfmethod(6) = 1;
    end
end 

for iarray_geo=1:length(FigConfig.arrays)
    if strcmp(FigConfig.arrays{iarray_geo}, 'all')
        do.array(:) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'linear_bs')
        do.array(1) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'linear_ef')
        do.array(2) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'rectangular')
        do.array(3) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'circular')    
        do.array(4) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'dualcircular')    
        do.array(5) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'stacked_circular')    
        do.array(6) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'spherical')
        do.array(7) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'circular_rigid_cylinder')    
        do.array(8) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'stacked_circular_rigid_cylinder')    
        do.array(9) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'circular_rigid_sphere')    
        do.array(10) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'spherical_rigid_sphere')    
        do.array(11) = 1; 
    end
end     

ArrayGeoInitials = {'LB','LE','R','C','DC','SC','S','C-RC','SC-RC','C-RS','S-RS'};
count=0;
for iarray_geo=1:narray_geo
    array_geo{iarray_geo} = ConfigAllArrays{iarray_geo,1,1}.Array.Geo;
    if do.array(iarray_geo)
        count=count+1;
        %array_geo_disp{count} = ConfigAllArrays{iarray_geo,1,1}.Array.GeoDisp;
        array_geo_disp{count} = ArrayGeoInitials{iarray_geo};
    end
end

%%% plotting beamforming measures in azimuth

if strcmp(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'3D') || strcmp(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'2Daz')


% Directivity index
if do.metric(1)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    %figure
                    %figure('pos',[10 10 1000 797]);
                    figure('pos',[5 5 700 550]);                  
                    bar(squeeze(BFEvalAllArrays.favg.DIaz(:,iM,ir,:,isource,inoise,ie))); 
                    %bar(DIazfavgAllArraysMethods{isource,inoise}(:,:,iM,ir)'); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Directivity Index in azimuth, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    
                    ylim([0 15]);
                    ylabel('DI_{az} [dB]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_DI_az_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

% White Noise Gain
if do.metric(2)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure('pos',[10 10 700 550]);
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.WNG(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['White Noise Gain in azimuth, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    
                    ylim([0 20]);
                    ylabel('WNG_{az} [dB]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_WNG_az_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end



% Beam width -3dB
if do.metric(3)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure('pos',[10 10 700 550]);
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.BWaz(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Beam Width -3dB in azimuth, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    
                    ylim([0 70]);
                    ylabel('BW_{-3dB az} [\circ]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_BW_az_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end


% Beam width -15dB
if do.metric(4)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure('pos',[10 10 700 550]);
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.BW15az(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Beam Width -15dB in azimuth, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    
                    %ylim([0 20]);
                    ylabel('BW_{-15dB az} [\circ]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_BW15_az_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

% Sidelobe suppresion level
if do.metric(5)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure('pos',[10 10 700 550]);
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.SLaz(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Sidelobe Suppresion Level in azimuth, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);                    
                    
                    ylim([0 15]);
                    ylabel('SSL_{az} [dB]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_SSL_az_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end


% Acoustic contrast
if do.metric(6)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure('pos',[10 10 700 550]);
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.Contrastaz(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Acoustic Contrast in azimuth, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);                    
                    
                    %ylim([0 20]);
                    ylabel('AC_{az} [\circ]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_AC_az_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

% Acoustic contrast Range
if do.metric(7) && isfield(BFEvalAllArrays,'ContrastRangeaz')
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure('pos',[10 10 700 550]);
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.ContrastRangeaz(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Acoustic Contrast Range in azimuth, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);                   
                    
                    %ylim([0 20]);
                    ylabel('AC Range_{az} [\circ]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_ACRange_az_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

end

if strcmp(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'3D') || strcmp(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'2Del')
    
% Directivity index
if do.metric(1)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.DIel(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Directivity Index in elevation, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    ylabel('DI_{el} [dB]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    %ylim([0 20]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_DI_el_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

% White Noise Gain
if do.metric(2)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.WNG(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['White Noise Gain in elevation, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    ylabel('WNG_{el} [dB]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    %ylim([0 20]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_WNG_el_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

% Beam width
if do.metric(3)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.BWel(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Beam Width in elevation, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    ylabel('BW_{el} [\circ]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    %ylim([0 20]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_BW_el_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

% Beam width -15dB
if do.metric(4)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.BW15el(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Beam Width -15dB in elevation, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    ylabel('BW_{-15dB el} [\circ]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    %ylim([0 20]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_BW15_el_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

% Sidelobe suppresion level
if do.metric(5)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.SLel(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Sidelobe Suppresion Level in elevation, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    ylabel('SSL_{el} [dB]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    %ylim([0 20]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_SSL_el_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

% Acoustic contrast
if do.metric(6)
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure
                    %figure('pos',[10 10 1000 797]);
                    bar(squeeze(BFEvalAllArrays.favg.Contrastel(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Acoustic Contrast in elevation, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    ylabel('AC_{el} [\circ]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    %ylim([0 20]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_AC_el_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

% Acoustic contrast
if do.metric(7) && isfield(BFEvalAllArrays,'ContrastRangeel')
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    %figure('pos',[10 10 1000 797]);
                    figure
                    bar(squeeze(BFEvalAllArrays.favg.ContrastRangeel(:,iM,ir,:,isource,inoise,ie))); 
                    set(gca,'XTickLabel',array_geo_disp);
                    %xticklabel_rotate([],FigConfig.Xtickrotangle,[],'fontsize',FigConfig.Fontsize);
                    %set(gca,'fontsize',FigConfig.Fontsize);
                    title(['Acoustic Contrast in elevation, ' ...
                            'M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    ylabel('AC_{el} [\circ]');
                    ylabelpos = get(get(gca,'Ylabel'),'Position');
                    set(get(gca,'Ylabel'),'Position',[ylabelpos(1)-0.05,ylabelpos(2),ylabelpos(3)]);
                    %ylim([0 20]);
                    legend(methodStr);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'Bar_AC_el_all_methods_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance '.eps']);
                    end
                end
            end
        end
    end
end
end

end
end
