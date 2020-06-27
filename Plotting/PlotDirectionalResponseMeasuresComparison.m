function PlotDirectionalResponseMeasuresComparison(ConfigAllArrays1,BFEvalAllArrays1, BFEvalAllArrays2, FigConfig)
% function PlotDirectionalResponseMeasuresComparison(ConfigAllArrays1,BFEvalAllArrays1, BFEvalAllArrays2, FigConfig)
% Function that plots a comparison of beamforming metrics between two sets 
% of simulations.
%
% input arguments:
%   ConfigAllArrays1: cell array containing all the configuration settings 
%       for one set of the simulations.
%   BFEvalAllArrays1: struct containing containing several performance
%       metrics for one set of the simulations.
%   BFEvalAllArrays2: struct containing containing several performance
%       metrics for the other set the simulations.
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.
%
% Note BFEvalAllArrays1 and BFEvalAllArrays2 are structs from the output of 
% rearrange_directional_response_measures.m. 


f = ConfigAllArrays1{1,1,1}.Filter.FreqVector;
blue = [0 0.447 0.741];
orange = [0.85,0.325,0.098];
red = [1 0 0];
yellow = [0.929 0.694 0.125];
purple = [0.494 0.184 0.556];
green = [0.4660 0.6740 0.188];
burdeaux = [0.635 0.447 0.741];
black = [0 0 0];

ArrayGeoInitials = {'LB','LE','R','C','DC','SC','S','C-RC','SC-RC','C-RS','S-RS'};
%for main geometries
%arrays with same geometry of microphones have same colour.
%arrays with the same baffle type have the same line style.
colours = {yellow,red,blue,purple,green,blue,blue,green};
lines = {'-','-','-','-','-','--','-.','-.'};

%for main geometries including stacked circular arrays
%arrays with same geometry of microphones have same colour.
%arrays with the same baffle type have the same line style.
colours = {yellow,orange,red,blue,purple,black,green,blue,black,blue,green};
lines = {'-','-','-','-','-','-','-','--','--','-.','-.'};

%for multicircular arrays
%arrays with same geometry of microphones have same colour (e.g. circular,
%dual circular, tricircular, rectangular, etc)
%arrays with the same baffle type have the same line style (e.g. equiangle
%arrangement, equispacing per circumference or uniform xy grid)
% colours = {blue,purple,'c',green,purple,'c',green,yellow,black,red,orange,'m'};
% lines = {'-','-','-','-','--','--','--','-.',':','-.','-.','-.'};

if isfield(BFEvalAllArrays1,'DIaz')
    plane = 'az';
elseif isfield(BFEvalAllArrays1,'DIel')
    plane = 'el';
end

%method to include part of a string as a field of a struct by using sprintf
%of the concatenated string ("DI"=part of variable + "plane"=newstring) and
%its output (string) include it in brackets to convert it to a field
narray_geo = size(BFEvalAllArrays1.(sprintf(strcat('DI',plane))),1);
nM = size(BFEvalAllArrays1.(sprintf(strcat('DI',plane))),2);
nr = size(BFEvalAllArrays1.(sprintf(strcat('DI',plane))),3);
nmethods = size(BFEvalAllArrays1.(sprintf(strcat('DI',plane))),4);
nsources = size(BFEvalAllArrays1.(sprintf(strcat('DI',plane))),5);
nnoises = size(BFEvalAllArrays1.(sprintf(strcat('DI',plane))),6);
nepsilons = size(BFEvalAllArrays1.(sprintf(strcat('DI',plane))),7);
nfreqs = length(f);

source1 = 1;
sourcen = 1;
noise1 = 1;
noisen = nnoises;
epsilon1 = 1;
epsilonn = 1;

M = zeros(1,nM);
r = zeros(1,nr);

for iM=1:nM
    M(iM) = size(ConfigAllArrays1{1,iM,1}.Array.MicPos,1);
end
%radius estimated from half the aperture of linear array
for ir=1:nr
    r(ir)= round(ConfigAllArrays1{1,1,ir}.Array.Aperture/2,2);
end

f=ConfigAllArrays1{1,1,1}.Filter.FreqVector;
foctaves = [32 63 125 250 500 1000 2000 4000 8000 16000];
foctaves = [63 125 250 500 1000 2000 4000 8000 16000];
ifirstoctave = find(foctaves>=f(1),1,'first');
ilastoctave = find(foctaves<=f(end),1,'last');
foctaves = foctaves(ifirstoctave:ilastoctave);
fbandcutoff = [foctaves.*2^(-1/2) foctaves(end).*2^(1/2)];
fLogPoints = length(foctaves)*2;
flogaxispos = evenLogAxis(f,fLogPoints);

if ~isfield(FigConfig,'NoiseType') || isempty(FigConfig.NoiseType)
    FigConfig.NoiseType=cell(1,nnoises);
else
    NoiseType=FigConfig.NoiseType;
    FigConfig.NoiseType=strcat('_',NoiseType);
end

if isfield(FigConfig,'Type') && strcmp(FigConfig.Type,'difference') 
    FigConfig.Error = 'error';
else
    FigConfig.Error = '';
end
    
rmm=r*1000;
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

array_geo = cell(1,narray_geo);

count=0;
for iarray_geo=1:narray_geo
    array_geo{iarray_geo} = ConfigAllArrays1{iarray_geo,1,1}.Array.Geo;
    if do.array(iarray_geo)
        count=count+1;
        %array_geo_disp{count} = ConfigAllArrays{iarray_geo,1,1}.Array.GeoDisp;
        array_geo_disp{count} = [ArrayGeoInitials{iarray_geo} ' ' FigConfig.Test1];
        count=count+1;
        array_geo_disp{count} = [ArrayGeoInitials{iarray_geo} ' ' FigConfig.Test2];
    end
end

%%% plotting beamforming measures in azimuth

if strcmp(ConfigAllArrays1{1,1,1}.ArrayMan.SteerSpace,'3D') || strcmp(ConfigAllArrays1{1,1,1}.ArrayMan.SteerSpace,'2Daz')

% Directivity index
if do.metric(1)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.DIaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.DIaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Directivity Index in azimuth, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('DI_{az} [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([0 20]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-10 10]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);
                    
                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'DI_az_' FigConfig.Error ConfigAllArrays1{1,1,1}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% White Noise Gain
if do.metric(2)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                        title(['White Noise Gain, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('WNG [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([-5 20]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-10 10]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'WNG_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end
 
% Beamwidth -3dB
if do.metric(3)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.BWaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.BWaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                        title(['Beam width in azimuth, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('BW_{az} [\circ]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([0 180]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-90 90]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'BW_az_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% Beamwidth -15dB
if do.metric(4)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.BW15az(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.BW15az(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Beam width -15dB in azimuth, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('BW_{-15dB az} [\circ]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([0 180]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-90 90]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'BW15_az_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% Sidelobe suppression level
if do.metric(5)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
   if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.SLaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.SLaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Sidelobe suppression level in azimuth, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('SSL_{az} [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([0 30]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-15 15]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'SSL_az_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
       end; end
end
end
end

% Acoustic Contrast
if do.metric(6)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.Contrastaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.Contrastaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Acoustic Contrast in azimuth, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('AC_{az} [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([0 50]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-25 25]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'AC_az_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% Acoustic Contrast within an angle range
if do.metric(7) && isfield(BFEvalAllArrays1,'ContrastRangeaz')
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.ContrastRangeaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.ContrastRangeaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Acoustic Contrast angle range in azimuth, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('AC_{range az} [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([-20 20]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-20 20]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'AC_range_az_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% On-axis Frequency Response
if do.metric(8) && isfield(BFEvalAllArrays1,'FRespaz')
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.FRespaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.FRespaz(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['On-axis Frequency Response in azimuth, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('20log_{10}|D(\phi_t,0)_{az}| [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([-30 10]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-15 15]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'FR_az_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

end


%%% plotting beamforming measures in elevation
if strcmp(ConfigAllArrays1{1,1,1}.ArrayMan.SteerSpace,'3D') || strcmp(ConfigAllArrays1{1,1,1}.ArrayMan.SteerSpace,'2Del')
    
% Directivity index
if do.metric(1)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.DIel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.DIel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Directivity Index in elevation, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('DI_{el} [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([0 20]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-10 10]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'DI_el_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% White Noise Gain
if do.metric(2)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.WNG(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['White Noise Gain, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('WNG [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([-5 20]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-10 10]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'WNG_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% Beamwidth
if do.metric(3)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.BWel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.BWel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Beam width in elevation, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('BW_{el} [\circ]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([0 180]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-90 90]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'BW_el_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% Beamwidth -15dB
if do.metric(4)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.BW15el(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.BW15el(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Beam width -15dB in elevation, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('BW_{-15dB el} [\circ]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([0 180]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-90 90]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'BW15_el_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% Sidelobe suppression level
if do.metric(5)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.SLel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.SLel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Sidelobe suppression level in elevation, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('SSL_{el} [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([0 30]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-15 15]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'SSL_el_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% Acoustic Contrast
if do.metric(6)
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.Contrastel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.Contrastel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Acoustic Contrast in elevation, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('AC_{el} [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([0 50]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-25 25]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'AC_el_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% Acoustic Contrast within an angle range
if do.metric(7) && isfield(BFEvalAllArrays1,'ContrastRangeel')
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.ContrastRangeel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.ContrastRangeel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['Acoustic Contrast angle range in elevation, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('AC_{range el} [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([-20 20]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-20 20]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'AC_range_el_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end

% On-axis Frequency Response
if do.metric(8) && isfield(BFEvalAllArrays1,'FRespel')
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    figure;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        semilogx(f,squeeze(BFEvalAllArrays1.FRespel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                        semilogx(f,squeeze(BFEvalAllArrays2.FRespel(iarray_geo,iM,ir,imethod,isource,inoise,ie,:)),...
                            'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo},'Marker','d','MarkerIndices',flogaxispos); hold on;
                        end
                    end
                    if FigConfig.Save; title(''); else
                    title(['On-axis Frequency Response in elevation, ' ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    xlabel('f [Hz]');
                    ylabel('20log_{10}|D(0,\theta_t)_{el}| [dB]');
                    if strcmp(FigConfig.Type,'absolute'); ylim([-30 10]);
                    elseif strcmp(FigConfig.Type,'difference'); ylim([-15 15]); end
                    legend(array_geo_disp);
                    %xticks(foctaves);
                    ax = gca;
                    ax.XTick = foctaves;
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);

                    if FigConfig.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            'FR_el_' FigConfig.Error ConfigAllArrays1{iarray_geo,iM,ir}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays1{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
end
end


end
end