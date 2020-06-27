function PlotBFWeights(ConfigAllArrays, BFWeightsAllArrays, FigConfig)
% function PlotBFWeights(ConfigAllArrays, BFWeightsAllArrays, FigConfig)
% function to plot the weights of different arrays and beamformers
% based on the settings in FigConfig. 
%
% input arguments:
%   ConfigAllArrays: cell array containing all the configuration settings 
%       for each microphone array, beamforming method, etc.
%   BFWeightsAllArrays: cell array containing several performance metrics 
%       for each array, beamforming method, etc. 
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.


FigConfig = defaultPlotValues(FigConfig);
[FigConfig,ArrayGeoInitials,BFInitials] = defaultArrayBFLineStyle(FigConfig);
if ~isfield(FigConfig,'FontSize') || isempty(FigConfig.FontSize); FigConfig.FontSize = 13; end

%method to include part of a string as a field of a struct by using sprintf
%of the concatenated string ("DI"=part of variable + "plane"=newstring) and
%its output (string) include it in brackets to convert it to a field
narray_geo = size(BFWeightsAllArrays,1);
nM = size(BFWeightsAllArrays,2);
nr = size(BFWeightsAllArrays,3);
nmethods = size(BFWeightsAllArrays{1,1,1},1);
nsources =  size(BFWeightsAllArrays{1,1,1},2);
nepsilons = 1;

source1 = 1;
sourcen = 1;
epsilon1 = 1;
epsilonn = 1;

f=ConfigAllArrays{1,1,1}.Filter.FreqVector;
t = 0:1/ConfigAllArrays{1,1,1}.Filter.Fs:(ConfigAllArrays{1,1,1}.Filter.Nfft-1)/ConfigAllArrays{1,1,1}.Filter.Fs;
samples = 1:ConfigAllArrays{1,1,1}.Filter.Nfft;
foctaves = [63 125 250 500 1000 2000 4000 8000 16000];
foctavesTicks = {'63', '125', '250', '500', '1k', '2k', '4k', '8k', '16k'};
ifirstoctave = find(foctaves>=f(1),1,'first');
ilastoctave = find(foctaves<=f(end),1,'last');
foctaves = foctaves(ifirstoctave:ilastoctave);
fbandcutoff = [foctaves.*2^(-1/2) foctaves(end).*2^(1/2)];


do.bfmethod = zeros(nmethods,1);
do.array = zeros(narray_geo,1);
if ~isfield(FigConfig,'WDomain')
    warning('Domain not specified. Plotting weights in frequency and time');
    FigConfig.WDomain = {'fd','td'};
end
do.domain = zeros(length(FigConfig.WDomain),1);
for idomain=1:length(FigConfig.WDomain)
    if strcmpi(FigConfig.WDomain{idomain}, 'fd'); do.domain(idomain) = true;
    elseif strcmpi(FigConfig.WDomain{idomain}, 'td'); do.domain(idomain) = true;
    else, error(['non-existent domain specified. Choose from ''','td','''or''','fd','''']);
    end
end


count = 1;
for imethod=1:length(FigConfig.Methods)
    if strcmpi(FigConfig.Methods{imethod}, 'all')
        do.bfmethod(:) = 1;
        methods = BFInitials;
    elseif strcmpi(FigConfig.Methods{imethod}, 'ds')
        %do.bfmethod(1) = 1;
        methods{count} = BFInitials{1};
    elseif strcmpi(FigConfig.Methods{imethod}, 'sda')
        %do.bfmethod(2) = 1;
        methods{count} = BFInitials{2};
    elseif strcmpi(FigConfig.Methods{imethod}, 'mvdr')
        %do.bfmethod(3) = 1;
        methods{count} = BFInitials{3};
    elseif strcmpi(FigConfig.Methods{imethod}, 'lcmv')    
        %do.bfmethod(4) = 1;
        methods{count} = BFInitials{5};
    elseif strcmpi(FigConfig.Methods{imethod}, 'ls')    
        %do.bfmethod(5) = 1;
        methods{count} = BFInitials{7};
    elseif strcmpi(FigConfig.Methods{imethod}, 'mvdr2')
        %do.bfmethod(6) = 1;
        methods{count} = BFInitials{4};
    elseif strcmpi(FigConfig.Methods{imethod}, 'lcmv2')
        %do.bfmethod(7) = 1;
        methods{count} = BFInitials{6};
    end
    
    count = count + 1;
    
    for iBF = 1:length(ConfigAllArrays{1,1,1}.Filter.Method)
        if strcmpi(ConfigAllArrays{1,1,1}.Filter.Method{iBF}, FigConfig.Methods{imethod})
            do.bfmethod(iBF) = 1;
        end
    end
end 

count=0;
if strcmpi(FigConfig.Arrays{1}, 'all')
    do.array(:) = 1; 
    ArrInd = 1:narray_geo;
else
    ArrInd = zeros(narray_geo,1);
    for iarray_geo=1:length(FigConfig.Arrays)
        for iArr = 1:narray_geo
        if strcmpi(ConfigAllArrays{iArr,1,1}.Array.GeoDisp, FigConfig.Arrays{iarray_geo})
            do.array(iArr) = 1;
            ArrInd(iarray_geo) = iArr;
        end
        end
    end
end
array_geo_disp = cell(length(FigConfig.Arrays),1);
count=0;
for iArr = 1:narray_geo; if do.array(iArr)
    count= count+1;
    for iallarray=1:length(ArrayGeoInitials)
        if strcmpi(ConfigAllArrays{iArr,1,1}.Array.GeoDisp,ArrayGeoInitials{iallarray})
            array_geo_disp{count} = ArrayGeoInitials{iallarray};
        end
    end
end; end     
 

%%% Plotting
iarr = 0;

for ie=epsilon1:epsilonn
if strcmpi(FigConfig.Compare,'arrays')
for imethod=1:nmethods; if do.bfmethod(imethod)
    for isource=source1:sourcen
        for iM=1:nM
            for ir=1:nr
                if FigConfig.Do.NewFig; figure; end
                for idomain=1:length(FigConfig.WDomain)
                    if strcmpi(FigConfig.WDomain{idomain},'fd')
                        for iarray_geo=1:narray_geo; if do.array(iarray_geo)  
                            w = squeeze(BFWeightsAllArrays{iarray_geo,iM,ir}{imethod,isource}.fd);
                            ax=subplot(length(FigConfig.WDomain)+1,1,1);
                            semilogx(f,db(w),'Color',FigConfig.LineStyle.Colours{iarray_geo},...
                                'LineStyle',FigConfig.LineStyle.Lines{iarray_geo}); hold on;
                            ax2=subplot(length(FigConfig.WDomain)+1,1,2);
                            semilogx(f,angle(w),'Color',FigConfig.LineStyle.Colours{iarray_geo},...
                                'LineStyle',FigConfig.LineStyle.Lines{iarray_geo}); hold on;
%                             ax=p1.ax;
%                             ax2=p2.ax;
                        end; end
                        if FigConfig.Do.Title
                            MethodStr = ConfigAllArrays{1,1,1}.Filter.Method{imethod};
                            if strcmpi(MethodStr,'ls'); MethodStr = [MethodStr ' N=' num2str(ConfigAllArrays{1,1,1}.Filter.ls.order(imethod))]; end
                            ax.Title.String = ['Weights in ' FigConfig.WDomain{idomain} ' for ' MethodStr];
                        end

                        ax.XLabel.Interpreter = 'latex';
                        ax.YLabel.Interpreter = 'latex';
                        ax.TickLabelInterpreter = 'latex';
                        if FigConfig.Do.XLabel; ax.XLabel.String = 'Frequency (Hz)'; end
                        if FigConfig.Do.YLabel; ax.YLabel.String = '$|w\left(\omega\right)|^2$ (dB)'; end
                        %ax.YLim = FigConfig.Temp.Ylim;
                        if FigConfig.Do.Legend; legend(ax,array_geo_disp); end
                        if FigConfig.Do.XTicks
                            ax.XTick = foctaves;
                            ax.XTickLabel = foctavesTicks;
                        end
                        ax.XLim = [max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))];
                        ax.XGrid = 'on'; ax.YGrid = 'on';
                        ax.FontSize = FigConfig.FontSize;
                        
                        ax2.XLabel.Interpreter = 'latex';
                        ax2.YLabel.Interpreter = 'latex';
                        ax2.TickLabelInterpreter = 'latex';
                        if FigConfig.Do.XLabel; ax2.XLabel.String = 'Frequency (Hz)'; end
                        if FigConfig.Do.YLabel; ax2.YLabel.String = '$\angle w\left(\omega\right)$ (rad)'; end
                        %ax.YLim = FigConfig.Temp.Ylim;
                        if FigConfig.Do.Legend; legend(ax2,array_geo_disp); end
                        if FigConfig.Do.XTicks
                            ax2.XTick = foctaves;
                            ax2.XTickLabel = foctavesTicks;
                        end
                        ax2.XLim = [max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))];
                        ax2.XGrid = 'on'; ax2.YGrid = 'on';
                        ax2.FontSize = FigConfig.FontSize;

                    elseif strcmpi(FigConfig.WDomain{idomain},'td')                        
                        for iarray_geo=1:narray_geo; if do.array(iarray_geo)  
                            w = squeeze(BFWeightsAllArrays{iarray_geo,iM,ir}{imethod,isource}.td);
                            ax=subplot(length(FigConfig.WDomain)+1,1,length(FigConfig.WDomain)*2-1);
                            plot(samples,w,'Color',FigConfig.LineStyle.Colours{iarray_geo},...
                                'LineStyle',FigConfig.LineStyle.Lines{iarray_geo}); hold on;
                        end; end
                        if FigConfig.Do.Title
                            MethodStr = ConfigAllArrays{1,1,1}.Filter.Method{imethod};
                            if strcmpi(MethodStr,'ls'); MethodStr = [MethodStr ' N=' num2str(ConfigAllArrays{1,1,1}.Filter.ls.order(imethod))]; end
                            ax.Title.String = ['Weights in ' FigConfig.WDomain{idomain} ' for ' MethodStr];
                        end

                        ax.XLabel.Interpreter = 'latex';
                        ax.YLabel.Interpreter = 'latex';
                        ax.TickLabelInterpreter = 'latex';
                        if FigConfig.Do.XLabel; ax.XLabel.String = 'Samples ($n$)'; end
                        if FigConfig.Do.YLabel; ax.YLabel.String = '$w\left(n\right)$'; end
                        %ax.YLim = FigConfig.Temp.Ylim;
                        if FigConfig.Do.Legend; legend(ax,array_geo_disp); end
                        ax.XLim = [samples(1) samples(end)];
                        ax.XGrid = 'on'; ax.YGrid = 'on';
                        ax.FontSize = FigConfig.FontSize;
                                               
                    end
                end   
                if FigConfig.Do.Save==1
                    print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                        'Weights_' ConfigAllArrays{1,1,1}.Filter.Method{imethod} 'all_arrays.eps']);
                end
            end
        end
    end
    end; end
elseif strcmpi(FigConfig.Compare,'methods')
iarr = 0;
for iarray_geo=1:narray_geo; if do.array(iarray_geo) 
    iarr = iarr + 1;
    for isource=source1:sourcen
        for iM=1:nM
            for ir=1:nr
                if FigConfig.Do.NewFig; figure; end
                for idomain=1:length(FigConfig.WDomain)
                    if strcmpi(FigConfig.WDomain{idomain},'fd')
                        for imethod=1:nmethods; if do.bfmethod(imethod)        
                            w = squeeze(BFWeightsAllArrays{iarray_geo,iM,ir}{imethod,isource}.fd);
                            ax=subplot(length(FigConfig.WDomain)+1,1,1);
                            semilogx(f,db(w),'Color',FigConfig.LineStyle.Colours{imethod},...
                                'LineStyle',FigConfig.LineStyle.Lines{imethod}); hold on;
                            ax2=subplot(length(FigConfig.WDomain)+1,1,2);
                            semilogx(f,angle(w),'Color',FigConfig.LineStyle.Colours{imethod},...
                                'LineStyle',FigConfig.LineStyle.Lines{imethod}); hold on;
%                             ax=p1.ax;
%                             ax2=p2.ax;
                        end; end
                        if FigConfig.Do.Title
                            ax.Title.String = ['Weights in ' FigConfig.WDomain{idomain} ' for ' ConfigAllArrays{iarray_geo,1,1}.Array.GeoDisp];
                        end

                        ax.XLabel.Interpreter = 'latex';
                        ax.YLabel.Interpreter = 'latex';
                        ax.TickLabelInterpreter = 'latex';
                        if FigConfig.Do.XLabel; ax.XLabel.String = 'Frequency (Hz)'; end
                        if FigConfig.Do.YLabel; ax.YLabel.String = '$|w\left(\omega\right)|^2$ (dB)'; end
                        %ax.YLim = FigConfig.Temp.Ylim;
                        if FigConfig.Do.Legend; legend(ax,methods); end
                        if FigConfig.Do.XTicks
                            ax.XTick = foctaves;
                            ax.XTickLabel = foctavesTicks;
                        end
                        ax.XLim = [max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))];
                        ax.XGrid = 'on'; ax.YGrid = 'on';
                        ax.FontSize = FigConfig.FontSize;
                        
                        ax2.XLabel.Interpreter = 'latex';
                        ax2.YLabel.Interpreter = 'latex';
                        ax2.TickLabelInterpreter = 'latex';
                        if FigConfig.Do.XLabel; ax2.XLabel.String = 'Frequency (Hz)'; end
                        if FigConfig.Do.YLabel; ax2.YLabel.String = '$\angle w\left(\omega\right)$ (rad)'; end
                        %ax.YLim = FigConfig.Temp.Ylim;
                        if FigConfig.Do.Legend; legend(ax2,methods); end
                        if FigConfig.Do.XTicks
                            ax2.XTick = foctaves;
                            ax2.XTickLabel = foctavesTicks;
                        end
                        ax2.XLim = [max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))];
                        ax2.XGrid = 'on'; ax2.YGrid = 'on';
                        ax2.FontSize = FigConfig.FontSize;

                    elseif strcmpi(FigConfig.WDomain{idomain},'td')                        
                        for imethod=1:nmethods; if do.bfmethod(imethod)                    
                            w = squeeze(BFWeightsAllArrays{iarray_geo,iM,ir}{imethod,isource}.td);
                            ax=subplot(length(FigConfig.WDomain)+1,1,length(FigConfig.WDomain)*2-1);
                            plot(samples,w,'Color',FigConfig.LineStyle.Colours{imethod},...
                                'LineStyle',FigConfig.LineStyle.Lines{imethod}); hold on;
                        end; end
                        if FigConfig.Do.Title
                            ax.Title.String = ['Weights in ' FigConfig.WDomain{idomain} ' for ' ConfigAllArrays{iarray_geo,1,1}.Array.GeoDisp];
                        end

                        ax.XLabel.Interpreter = 'latex';
                        ax.YLabel.Interpreter = 'latex';
                        ax.TickLabelInterpreter = 'latex';
                        if FigConfig.Do.XLabel; ax.XLabel.String = 'Samples ($n$)'; end
                        if FigConfig.Do.YLabel; ax.YLabel.String = '$w\left(n\right)$'; end
                        %ax.YLim = FigConfig.Temp.Ylim;
                        if FigConfig.Do.Legend; legend(ax2,methods); end
                        ax.XLim = [samples(1) samples(end)];
                        ax.XGrid = 'on'; ax.YGrid = 'on';
                        ax.FontSize = FigConfig.FontSize;
                                               
                    end
                end

                if FigConfig.Do.Save==1
                    print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                        'Weights_all_methods_' ConfigAllArrays{iarray_geo,1,1}.Array.GeoDisp '.eps']);
                end

            end
        end
    end
end; end
end
end

end

