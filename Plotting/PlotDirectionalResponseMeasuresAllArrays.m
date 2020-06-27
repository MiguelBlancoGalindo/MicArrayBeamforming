function [FigConfig,do] = PlotDirectionalResponseMeasuresAllArrays(ConfigAllArrays,BFEvalAllArrays, FigConfig)
% function [FigConfig,do] = PlotDirectionalResponseMeasuresAllArrays(ConfigAllArrays,BFEvalAllArrays, FigConfig)
% Function that plots performance metrics from the beamforming analysis.
%
% input arguments:
%   ConfigAllArrays: cell array containing all the configuration settings 
%       for each microphone array, beamforming method, etc.
%   BFEvalAllArrays: struct containing several performance metrics for each 
%       array, beamforming method, etc. Note BFEvalAllArrays is a struct 
%       from the output of rearrange_directional_response_measures.m. 
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.
%
% output arguments:
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.
%   do: struct contaning logic statements of which metrics, arrays and
%       beamformers to plot.



if ~isfield(FigConfig,'FontSize') || isempty(FigConfig.FontSize); FigConfig.FontSize = 13; end
f = ConfigAllArrays{1,1,1}.Filter.FreqVector;
blue = [0 0.447 0.741];
orange = [0.85,0.325,0.098];
red = [1 0 0];
yellow = [0.929 0.694 0.125];
purple = [0.494 0.184 0.556];
green = [0.4660 0.6740 0.188];
burdeaux = [0.635 0.447 0.741];
black = [0 0 0];
darkGreen = [0.21 0.69 0.04];
ArrayGeoInitials = {'LB','LE','R','C','DC','SC','S','C-RC','SC-RC','C-RS','S-RS'};
BFInitials = {'DSB','SDB','MVDRB1','MVDRB2','LCMVB1','LCMVB2','LSB'};

if strcmpi(FigConfig.Compare,'arrays')
%for main geometries including stacked circular arrays
%arrays with same geometry of microphones have same colour.
%arrays with the same baffle type have the same line style.
FigConfig.Colours = {yellow,orange,red,blue,purple,black,green,blue,black,blue,green};
FigConfig.Lines = {'-','-','-','-','-','-','-','--','--','-.','-.'};
%for multicircular arrays
%arrays with same geometry of microphones have same colour (e.g. circular,
%dual circular, tricircular, rectangular, etc)
%arrays with the same baffle type have the same line style (e.g. equiangle
%arrangement, equispacing per circumference or uniform xy grid)
% colours = {blue,purple,'c',green,purple,'c',green,yellow,black,red,orange,'m'};
% lines = {'-','-','-','-','--','--','--','-.',':','-.','-.','-.'};

elseif strcmpi(FigConfig.Compare,'methods')  
    FigConfig.Colours = {yellow,orange,red,blue,purple,black,green,burdeaux,darkGreen};
    for icolour = 1:length(FigConfig.Colours)
        FigConfig.Lines{icolour} = '-';
        FigConfig.Markers{icolour} = 'none';
    end
end
   
if isfield(BFEvalAllArrays,'DIaz')
    plane = 'az';
elseif isfield(BFEvalAllArrays,'DIel')
    plane = 'el';
elseif isfield(BFEvalAllArrays,'DI3d')
    plane = '3d';
end
FigConfig.SteerSpace = ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace;
%method to include part of a string as a field of a struct by using sprintf
%of the concatenated string ("DI"=part of variable + "plane"=newstring) and
%its output (string) include it in brackets to convert it to a field
narray_geo = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),1);
nM = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),2);
nr = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),3);
nmethods = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),4);
nsources = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),5);
nnoises = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),6);
nepsilons = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),7);
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
    M(iM) = size(ConfigAllArrays{1,iM,1}.Array.MicPos,1);
end
%radius estimated from half the aperture of linear array
for ir=1:nr
    r(ir)= round(ConfigAllArrays{1,1,ir}.Array.Aperture/2,2);
end

f=ConfigAllArrays{1,1,1}.Filter.FreqVector;
foctaves = [63 125 250 500 1000 2000 4000 8000 16000];
foctavesTicks = {'63', '125', '250', '500', '1k', '2k', '4k', '8k', '16k'};
ifirstoctave = find(foctaves>=f(1),1,'first');
ilastoctave = find(foctaves<=f(end),1,'last');
foctaves = foctaves(ifirstoctave:ilastoctave);
fbandcutoff = [foctaves.*2^(-1/2) foctaves(end).*2^(1/2)];

if ~isfield(FigConfig,'NoiseType') || isempty(FigConfig.NoiseType)
    FigConfig.NoiseType=cell(1,nnoises);
else
    NoiseType=FigConfig.NoiseType;
    FigConfig.NoiseType=strcat('_',NoiseType);
end

if isfield(FigConfig,'Type') && strcmpi(FigConfig.Type,'difference') 
    FigConfig.Error = 'error';
else
    FigConfig.Error = '';
end
if ~isfield(FigConfig,'Units')
    FigConfig.Units = 'normalized';
end
rmm=r*1000;
nmetrics = 9; 
ntotalmethods = nmethods;
do.metric = zeros(nmetrics,1);
do.bfmethod = zeros(ntotalmethods,1);
do.array = zeros(narray_geo,1);
metricStr = cell(nmetrics,1);
for imetric=1:length(FigConfig.Metrics)
    if strcmpi(FigConfig.Metrics{imetric}, 'all')
        do.metric(:) = 1;
    elseif strcmpi(FigConfig.Metrics{imetric}, 'BW')
        do.metric(1) = 1;
        metricStr{1} = (sprintf(strcat(FigConfig.Metrics{imetric},plane)));
    elseif strcmpi(FigConfig.Metrics{imetric}, 'BW15')
        do.metric(2) = 1;
        metricStr{2} = (sprintf(strcat(FigConfig.Metrics{imetric},plane)));
    elseif strcmpi(FigConfig.Metrics{imetric}, 'SSL')    
        do.metric(3) = 1;
        metricStr{3} = (sprintf(strcat('SL',plane)));
    elseif strcmpi(FigConfig.Metrics{imetric}, 'AC')    
        do.metric(4) = 1;
        metricStr{4} = (sprintf(strcat('Contrast',plane)));
    elseif strcmpi(FigConfig.Metrics{imetric}, 'ACRange')
        do.metric(5) = 1;
        metricStr{5} = (sprintf(strcat('ContrastRange',plane)));
    elseif strcmpi(FigConfig.Metrics{imetric}, 'DI')
        do.metric(6) = 1;
        metricStr{6} = (sprintf(strcat(FigConfig.Metrics{imetric},plane)));
    elseif strcmpi(FigConfig.Metrics{imetric}, 'WNG')
        do.metric(7) = 1;
        metricStr{7} = 'WNG';
    elseif strcmpi(FigConfig.Metrics{imetric}, 'FRange')    
        do.metric(8) = 1;
        metricStr{8} = (sprintf(strcat(FigConfig.Metrics{imetric},plane)));
%     elseif strcmpi(FigConfig.Metrics{imetric}, 'BW3')
%         do.metric(9) = 1;
    end
    %metricStr{imetric} = (sprintf(strcat(FigConfig.Metrics{imetric},plane)));
end 

count = 1;
for imethod=1:length(FigConfig.Methods)
    if strcmpi(FigConfig.Methods{imethod}, 'all')
        do.bfmethod(:) = 1;
        methods = BFInitials;
    elseif strcmpi(FigConfig.Methods{imethod}, 'ds')
        methods{count} = BFInitials{1};
        if strcmpi(FigConfig.Compare,'methods')    
            FigConfig.Colours{count} = 'c';
            FigConfig.Lines{count} = '-';
            FigConfig.Markers{count} = 'none';
        end
    elseif strcmpi(FigConfig.Methods{imethod}, 'sda')
        methods{count} = BFInitials{2};
        if strcmpi(FigConfig.Compare,'methods')    
            FigConfig.Colours{count} = 'b';
            FigConfig.Lines{count} = '-';
            FigConfig.Markers{count} = 'none';
        end
    elseif strcmpi(FigConfig.Methods{imethod}, 'mvdr')
        methods{count} = BFInitials{3};
        if strcmpi(FigConfig.Compare,'methods')    
            FigConfig.Colours{count} = 'r';
            FigConfig.Lines{count} = '-.';
            FigConfig.Markers{count} = 'none';
        end
    elseif strcmpi(FigConfig.Methods{imethod}, 'lcmv')    
        methods{count} = BFInitials{5};
        if strcmpi(FigConfig.Compare,'methods')    
            FigConfig.Colours{count} = burdeaux;
            FigConfig.Lines{count} = '-';
            FigConfig.Markers{count} = 'none';
        end
    elseif strcmpi(FigConfig.Methods{imethod}, 'ls')    
        methods{count} = BFInitials{7};
        if strcmpi(FigConfig.Compare,'methods')    
            FigConfig.Colours{count} = darkGreen;
            FigConfig.Lines{count} = '-.';
            FigConfig.Markers{count} = 'none';
        end
    elseif strcmpi(FigConfig.Methods{imethod}, 'mvdr2')
        methods{count} = BFInitials{4};
        if strcmpi(FigConfig.Compare,'methods')    
            FigConfig.Colours{count} = 'r';
            FigConfig.Lines{count} = '--';
            FigConfig.Markers{count} = 'none';
        end
    elseif strcmpi(FigConfig.Methods{imethod}, 'lcmv2')
        methods{count} = BFInitials{6};
        if strcmpi(FigConfig.Compare,'methods')    
            FigConfig.Colours{count} = burdeaux;
            FigConfig.Lines{count} = '--';
            FigConfig.Markers{count} = 'none';
        end
    end
    
    count = count + 1;
    
    for iBF = 1:length(ConfigAllArrays{1,1,1}.Filter.Method)
        if strcmpi(ConfigAllArrays{1,1,1}.Filter.Method{iBF}, FigConfig.Methods{imethod})
            do.bfmethod(iBF) = 1;
        end
    end
end 

for iarray_geo=1:length(FigConfig.Arrays)
    if strcmpi(FigConfig.Arrays{iarray_geo}, 'all')
        do.array(:) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'linear_bs')
        do.array(1) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'linear_ef')
        do.array(2) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'rectangular')
        do.array(3) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'circular')    
        do.array(4) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'dualcircular')    
        do.array(5) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'stacked_circular')    
        do.array(6) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'spherical')
        do.array(7) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'circular_rigid_cylinder')    
        do.array(8) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'stacked_circular_rigid_cylinder')    
        do.array(9) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'circular_rigid_sphere')    
        do.array(10) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray_geo}, 'spherical_rigid_sphere')    
        do.array(11) = 1; 
    end
end     

array_geo = cell(1,narray_geo);

count=0;
for iarray_geo=1:narray_geo
    array_geo{iarray_geo} = ConfigAllArrays{iarray_geo,1,1}.Array.Geo;
    if do.array(iarray_geo)
        count=count+1;
        array_geo_disp{count} = ConfigAllArrays{iarray_geo,1,1}.Array.GeoDisp;
    end
end

%finding frequency range to display for each bfmethod and array according
%to FigConfig
%fRange = zeros(narray_geo,nmethods,2);
ifRange = zeros(narray_geo,nM,nr,nmethods,noisen,sourcen,epsilonn,2);
if ~isfield(FigConfig,'FRange') || isempty(FigConfig.FRange) || strcmpi(FigConfig.FRange,'fullband')
    ifRange = repmat(shiftdim([1; length(f)],-ndims(ifRange)),narray_geo,nM,nr,nmethods,noisen,sourcen,epsilonn,1);
    %fRange(iarray,imethod,:) = repmat(shiftdim(f,-1),narray_geo,nmethods,1);
else
    for ie=epsilon1:epsilonn
        for imethod=1:nmethods
            for isource=source1:sourcen
                for inoise=noise1:noisen
                    for iM=1:nM
                        for ir=1:nr
                            for iarray=1:narray_geo
                                if strcmpi(FigConfig.FRange,'operative')
                                    [~,ifRange(iarray,iM,ir,imethod,isource,inoise,ie,1)] = min(abs(f-squeeze(BFEvalAllArrays.(sprintf(strcat('FRange',plane)))(iarray,iM,ir,imethod,isource,inoise,ie,1))));
                                    [~,ifRange(iarray,iM,ir,imethod,isource,inoise,ie,2)] = min(abs(f-squeeze(BFEvalAllArrays.(sprintf(strcat('FRange',plane)))(iarray,iM,ir,imethod,isource,inoise,ie,2))));
                                    %fRange = f(ifRange(iarray,imethod,1):ifRange(iarray,imethod,2));
                                elseif strcmpi(FigConfig.FRange,'nonaliased')
                                    ifRange(iarray,iM,ir,imethod,isource,inoise,ie,1) = 1;
                                    [~,ifRange(iarray,iM,ir,imethod,isource,inoise,ie,2)] = min(abs(f-squeeze(BFEvalAllArrays.(sprintf(strcat('FRange',plane)))(iarray,iM,ir,imethod,isource,inoise,ie,2))));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


%%% Plotting
imet = 0;
for imetric=1:nmetrics; if do.metric(imetric)
imet = imet +1;
for ie=epsilon1:epsilonn
if strcmpi(FigConfig.Compare,'arrays')
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    if imet==1 || FigConfig.Do.NewFig; figure; end
                    if isfield(FigConfig.Do,'Subplot') && FigConfig.Do.Subplot
                        axes('Units',FigConfig.Units,'Position',[FigConfig.Subplot.Pos.Left(imet),FigConfig.Subplot.Pos.Bottom(imet),FigConfig.Subplot.Width,FigConfig.Subplot.Height]);
                    end
                    iarr = 0;
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                            iarr = iarr+1;
                            if1 = ifRange(iarray_geo,iM,ir,imethod,isource,inoise,ie,1);
                            if2 = ifRange(iarray_geo,iM,ir,imethod,isource,inoise,ie,2);
                            fRange = f(if1:if2);
                            metric = squeeze(BFEvalAllArrays.(sprintf(metricStr{imetric}))(iarray_geo,iM,ir,imethod,isource,inoise,ie,if1:if2));
                            if iarr==1 
                                p=semilogx(fRange,metric,...
                                'Color',FigConfig.Colours{iarray_geo},'LineStyle',FigConfig.Lines{iarray_geo}); hold on;
                                ax = p.Parent;
                            else
                                semilogx(ax,fRange,metric,...
                                    'Color',FigConfig.Colours{iarray_geo},'LineStyle',FigConfig.Lines{iarray_geo}); hold on;
                            end
                            
                            if isfield(FigConfig,'FRange') && ~strcmpi(FigConfig.FRange,'fullband')
                                if ~isfield(FigConfig,'Ylim')
                                    FigConfig.Ylim = [ax.YLim(1) ax.YLim(2)*1.2];
                                end
%                                 semilogx(ax,[fRange(1) fRange(1)],[FigConfig.Ylim(2)*0.9 FigConfig.Ylim(2)],...
%                                     'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
%                                 semilogx(ax,[fRange(end) fRange(end)],[FigConfig.Ylim(2)*0.9 FigConfig.Ylim(2)],...
%                                     'Color',colours{iarray_geo},'LineStyle',lines{iarray_geo}); hold on;
                            end

                            if isfield(FigConfig.Do,'Zoom') && FigConfig.Do.Zoom
                                if iarr==1
                                    ax2 = axes('Position',[FigConfig.Zoom.Pos.Left(imet),FigConfig.Zoom.Pos.Bottom(imet),FigConfig.Zoom.Width,FigConfig.Zoom.Height]);
                                end
                                semilogx(ax2,fRange,metric,...
                                    'Color',FigConfig.Colours{iarray_geo},'LineStyle',FigConfig.Lines{iarray_geo}); hold on;
                            end

                    end; end

                    FigConfig = getFigProperties(FigConfig, imetric);
                    if FigConfig.Do.Title
                        title([FigConfig.Temp.Title ' in ' plane ', ' ConfigAllArrays{iarray_geo,iM,ir}.Filter.Method{1,imethod} ...
                             ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end
                    ax.XLabel.Interpreter = 'latex';
                    ax.YLabel.Interpreter = 'latex';
                    ax.TickLabelInterpreter = 'latex';
                    if FigConfig.Do.XLabel; ax.XLabel.String = 'Frequency (Hz)'; end
                    if FigConfig.Do.YLabel; ax.YLabel.String = FigConfig.Temp.Ylabel; end
                    if ~isfield(FigConfig,'Ylim') || isempty(FigConfig.Ylim)
                        ax.YLim = FigConfig.Temp.Ylim;
                    end
                    if FigConfig.Do.Legend; legend(array_geo_disp); end
                    if FigConfig.Do.XTicks
                        ax.XTick = foctaves;
                        ax.XTickLabel = foctavesTicks;
                    end
                    ax.XLim = [max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))];
                    ax.XGrid = 'on'; ax.YGrid = 'on';
                    ax.FontSize = FigConfig.FontSize;
                    
                    if isfield(FigConfig.Do,'Zoom') && FigConfig.Do.Zoom
                        ax2.TickLabelInterpreter = 'latex';
                        ax2.YLim = FigConfig.Zoom.YLim;
                        if FigConfig.Do.XTicks
                            ax2.XTick = foctaves;
                            ax2.XTickLabel = foctavesTicks;
                        end
                        ax2.XLim = [FigConfig.Zoom.FRange(1) FigConfig.Zoom.FRange(2)];
                        ax2.XGrid = 'on'; ax2.YGrid = 'on';
                        ax2.FontSize = FigConfig.FontSize;
                    end
                    
                    if FigConfig.Do.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            MetricStr '_' plane '_' FigConfig.Error ConfigAllArrays{1,1,1}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
        end; end
end
elseif strcmpi(FigConfig.Compare,'methods')
for iarray_geo=1:narray_geo; if do.array(iarray_geo) 
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    if imet==1 || FigConfig.Do.NewFig; figure; end
                    if isfield(FigConfig.Do,'Subplot') && FigConfig.Do.Subplot 
                        ax=axes('Units',FigConfig.Units,'Position',[FigConfig.Subplot.Pos.Left(imet),FigConfig.Subplot.Pos.Bottom(imet),FigConfig.Subplot.Width,FigConfig.Subplot.Height]);
                    end
                    imeth = 0;
                    for imethod=1:nmethods; if do.bfmethod(imethod) 
                        imeth = imeth + 1;
                        if1 = ifRange(iarray_geo,iM,ir,imethod,isource,inoise,ie,1);
                        if2 = ifRange(iarray_geo,iM,ir,imethod,isource,inoise,ie,2);
                        fRange = f(if1:if2);
                        metric = squeeze(BFEvalAllArrays.(sprintf(metricStr{imetric}))(iarray_geo,iM,ir,imethod,isource,inoise,ie,if1:if2));
                        if imeth==1 
                            p=semilogx(fRange,metric,...
                            'Color',FigConfig.Colours{imethod},'LineStyle',FigConfig.Lines{imethod},'Marker',FigConfig.Markers{imethod}); hold on;
                            ax = p.Parent;
                        else
                            semilogx(ax,fRange,metric,...
                                'Color',FigConfig.Colours{imethod},'LineStyle',FigConfig.Lines{imethod},'Marker',FigConfig.Markers{imethod}); hold on;
                        end
                        
                        if isfield(FigConfig,'FRange') && ~strcmpi(FigConfig.FRange,'fullband')
                            if ~isfield(FigConfig,'Ylim')
                                FigConfig.Ylim = [ax.YLim(1) ax.YLim(2)*1.2];
                            end
%                             semilogx(ax,[fRange(1) fRange(1)],[FigConfig.Ylim(2)*0.9 FigConfig.Ylim(2)],...
%                                 'Color',colours{imethod},'LineStyle',lines{imethod}); hold on;
%                             semilogx(ax,[fRange(end) fRange(end)],[FigConfig.Ylim(2)*0.9 FigConfig.Ylim(2)],...
%                                 'Color',colours{imethod},'LineStyle',lines{imethod}); hold on;
                        end

                        if isfield(FigConfig.Do,'Zoom') && FigConfig.Do.Zoom
                            if imeth==1
                                ax2 = axes('Position',[FigConfig.Zoom.Pos.Left(imet),FigConfig.Zoom.Pos.Bottom(imet),FigConfig.Zoom.Width,FigConfig.Zoom.Height]);
                            end
                            semilogx(ax2,fRange,metric,...
                                'Color',FigConfig.Colours{imethod},'LineStyle',FigConfig.Lines{imethod}); hold on;
                        end
                        
                    end; end
                    
                    FigConfig = getFigProperties(FigConfig, imetric);
                    if FigConfig.Do.Title
                        title([FigConfig.Temp.Title ' in ' plane ', '  array_geo_disp{imeth} ...
                            ', M=' num2str(M(iM)) ', r=' num2str(r(ir)) 'm']);
                    end

                    ax.XLabel.Interpreter = 'latex';
                    ax.YLabel.Interpreter = 'latex';
                    ax.TickLabelInterpreter = 'latex';
                    if FigConfig.Do.XLabel; ax.XLabel.String = 'Frequency (Hz)'; end
                    if FigConfig.Do.YLabel; ax.YLabel.String = FigConfig.Temp.Ylabel; end
                    if ~isfield(FigConfig,'Ylim') || isempty(FigConfig.Ylim)
                        ax.YLim = FigConfig.Temp.Ylim;
                    end
                    if FigConfig.Do.Legend; legend(methods); end
                    if FigConfig.Do.XTicks
                        ax.XTick = foctaves;
                        ax.XTickLabel = foctavesTicks;
                    end
                    ax.XLim = [max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))];
                    ax.XGrid = 'on'; ax.YGrid = 'on';
                    ax.FontSize = FigConfig.FontSize;
                    
                    if isfield(FigConfig.Do,'Zoom') && FigConfig.Do.Zoom
                        ax2.TickLabelInterpreter = 'latex';
                        ax2.YLim = FigConfig.Zoom.YLim;
                        if FigConfig.Do.XTicks
                            ax2.XTick = foctaves;
                            ax2.XTickLabel = foctavesTicks;
                        end
                        ax2.XLim = [FigConfig.Zoom.FRange(1) FigConfig.Zoom.FRange(2)];
                        ax2.XGrid = 'on'; ax2.YGrid = 'on';
                        ax2.FontSize = FigConfig.FontSize;
                    end
                    
                    if FigConfig.Do.Save==1
                        print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                            MetricStr '_' plane '_' FigConfig.Error ConfigAllArrays{1,1,1}.Filter.Method{1,imethod} '_all_arrays_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{1,1,1}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                    end
                end
            end
        end
    end
end; end
end
end
end; end


end


function [FigConfig] = getFigProperties(FigConfig, imetric)
    if strcmpi(FigConfig.SteerSpace,'2Daz')
        SpaceVar = '\varphi';
    elseif strcmpi(FigConfig.SteerSpace,'2Del')
        SpaceVar = '\vartheta';
    elseif strcmpi(FigConfig.SteerSpace,'3D')
        SpaceVar = '\Phi';
    end
    if imetric == 1
        FigConfig.Temp.Title = 'Beam Width'; 
        FigConfig.Temp.Ylabel = ['BW($' SpaceVar '$) \ ($^\circ$)'];
        if strcmpi(FigConfig.Type,'absolute'); FigConfig.Temp.Ylim = [0 180];
        elseif strcmpi(FigConfig.Type,'difference'); FigConfig.Temp.Ylim = [-90 90]; end
    elseif imetric == 2
        FigConfig.Temp.Title = 'Beam Width -15dB';  
        FigConfig.Temp.Ylabel = ['BW($' SpaceVar '$) \ ($^\circ$)'];
        if strcmpi(FigConfig.Type,'absolute'); FigConfig.Temp.Ylim = [0 180];
        elseif strcmpi(FigConfig.Type,'difference'); FigConfig.Temp.Ylim = [-90 90]; end
    elseif imetric == 3
        FigConfig.Temp.Title = 'Sidelobe Suppression';   
        FigConfig.Temp.Ylabel = ['SSL($' SpaceVar '$) \ (dB)'];
        if strcmpi(FigConfig.Type,'absolute'); FigConfig.Temp.Ylim = [0 20];
        elseif strcmpi(FigConfig.Type,'difference'); FigConfig.Temp.Ylim = [-10 10]; end
    elseif imetric == 4
        FigConfig.Temp.Title = 'Acoustic Contrast';     
        FigConfig.Temp.Ylabel = ['AC($' SpaceVar '$) \ (dB)'];
        if strcmpi(FigConfig.Type,'absolute'); FigConfig.Temp.Ylim = [0 50];
        elseif strcmpi(FigConfig.Type,'difference'); FigConfig.Temp.Ylim = [-25 25]; end
    elseif imetric == 5
        FigConfig.Temp.Title = 'Acoustic Contrast Range';    
        FigConfig.Temp.Ylabel = ['AC$_{Range}$($' SpaceVar '$) \ (dB)'];
        if strcmpi(FigConfig.Type,'absolute'); FigConfig.Temp.Ylim = [0 20];
        elseif strcmpi(FigConfig.Type,'difference'); FigConfig.Temp.Ylim = [-20 20]; end
    elseif imetric == 6
        FigConfig.Temp.Title = 'Directivity Index';
        FigConfig.Temp.Ylabel = ['DI($' SpaceVar '$) \ (dB)'];
        if strcmpi(FigConfig.Type,'absolute'); FigConfig.Temp.Ylim = [0 20];
        elseif strcmpi(FigConfig.Type,'difference'); FigConfig.Temp.Ylim = [-10 10]; end
    elseif imetric == 7
        FigConfig.Temp.Title = 'White Noise Gain';
        FigConfig.Temp.Ylabel = ['WNG($' SpaceVar '$) \ (dB)'];
        if strcmpi(FigConfig.Type,'absolute'); FigConfig.Temp.Ylim = [-15 20];
        elseif strcmpi(FigConfig.Type,'difference'); FigConfig.Temp.Ylim = [-10 10]; end
    elseif imetric == 8
        FigConfig.Temp.Title = 'On-axis Frequency Response';       
        FigConfig.Temp.Ylabel = ['20log$_{10}$|d($' SpaceVar '_t$)| \ (dB)'];
        if strcmpi(FigConfig.Type,'absolute'); FigConfig.Temp.Ylim = [-30 10];
        elseif strcmpi(FigConfig.Type,'difference'); FigConfig.Temp.Ylim = [-15 15]; end
    end
end
