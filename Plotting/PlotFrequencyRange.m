function p = PlotFrequencyRange(BFEvalAllArrays,ConfigAllArrays, FRangeLS, FigConfig)
% function p = PlotFrequencyRange(BFEvalAllArrays,ConfigAllArrays, FRangeLS, FigConfig)
% function to plot the frequency range of different arrays and beamformers
% based on the settings in FigConfig. It also includes the frequency
% invariant range of the least-squares beamformer (LSB) delimited by 
% vertical lines.
%
% input arguments:
%   BFEvalAllArrays: struct containing several performance metrics for each 
%       array, beamforming method, etc. Note BFEvalAllArrays is a struct 
%       from the output of rearrange_directional_response_measures.m. 
%   ConfigAllArrays: cell array containing all the configuration settings 
%       for each microphone array, beamforming method, etc.
%   FRangeLS: frequency invariant range of the LSB. This can be obtained by
%       running directionalResponseErrorLSB.m.
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.
%
% output arguments:
%   p: plot handle.

if ~isfield(FigConfig,'FontSize') || isempty(FigConfig.FontSize); FigConfig.FontSize = 13; end
f = ConfigAllArrays{1,1,1}.Filter.FreqVector;

burdeaux = [0.635 0.447 0.741];
darkGreen = [0.21 0.69 0.04];
ArrayGeoInitials = {'LB','LE','R','C','DC','SC','S','C-RC','SC-RC','C-RS','S-RS'};
BFInitials = {'DSB','SDB','MVDRB1','MVDRB2','LCMVB1','LCMVB2','LSB'};
   
colours = {'c','b','r','r',burdeaux,burdeaux,darkGreen};
lines = {'-','-','-.','--','-','--','-.'};
markers = {'none','none','none','none','none','none','none'};
lineWidth = 1;

if isfield(BFEvalAllArrays,'DIaz')
    plane = 'az';
elseif isfield(BFEvalAllArrays,'DIel')
    plane = 'el';
end

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
foctaves = [32 63 125 250 500 1000 2000 4000 8000 16000];
foctavesTicks = {'32','63', '125', '250', '500', '1k', '2k', '4k', '8k', '16k'};
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
    
rmm=r*1000;
nmetrics = 9; 
ntotalmethods = nmethods;
do.BFMethod = zeros(ntotalmethods,1);
do.Array = zeros(narray_geo,1);

count = 1;
for imethod=1:length(FigConfig.Methods)
    if strcmpi(FigConfig.Methods{imethod}, 'all')
        do.BFMethod(:) = 1;
        methods = BFInitials;
    elseif strcmpi(FigConfig.Methods{imethod}, 'ds')
        %do.BFMethod(1) = 1;
        methods{count} = BFInitials{1};
        colours{count} = 'c';
        lines{count} = '-';
        markers{count} = 'none';
    elseif strcmpi(FigConfig.Methods{imethod}, 'sda')
        %do.BFMethod(2) = 1;
        methods{count} = BFInitials{2};
        colours{count} = 'b';
        lines{count} = '-';
        markers{count} = 'none';
    elseif strcmpi(FigConfig.Methods{imethod}, 'mvdr')
        %do.BFMethod(3) = 1;
        methods{count} = BFInitials{3};
        colours{count} = 'r';
        lines{count} = '-.';
        markers{count} = 'none';
    elseif strcmpi(FigConfig.Methods{imethod}, 'lcmv')    
        %do.BFMethod(4) = 1;
        methods{count} = BFInitials{5};
        colours{count} = burdeaux;
        lines{count} = '-';
        markers{count} = 'none';
    elseif strcmpi(FigConfig.Methods{imethod}, 'ls')    
        %do.BFMethod(5) = 1;
        methods{count} = BFInitials{7};
        colours{count} = darkGreen;
        lines{count} = '-.';
        markers{count} = 'none';
    elseif strcmpi(FigConfig.Methods{imethod}, 'mvdr2')
        %do.BFMethod(6) = 1;
        methods{count} = BFInitials{4};
        colours{count} = 'r';
        lines{count} = '--';
        markers{count} = 'none';
    elseif strcmpi(FigConfig.Methods{imethod}, 'lcmv2')
        %do.BFMethod(7) = 1;
        methods{count} = BFInitials{6};
        colours{count} = burdeaux;
        lines{count} = '--';
        markers{count} = 'none';
    end
    
    count = count + 1;
    
    for iBF = 1:length(ConfigAllArrays{1,1,1}.Filter.Method)
        if strcmpi(ConfigAllArrays{1,1,1}.Filter.Method{iBF}, FigConfig.Methods{imethod})
            do.BFMethod(iBF) = 1;
        end
    end
end 

for iarray=1:length(FigConfig.Arrays)
    if strcmpi(FigConfig.Arrays{iarray}, 'all')
        do.Array(:) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'linear_bs')
        do.Array(1) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'linear_ef')
        do.Array(2) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'rectangular')
        do.Array(3) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'circular')    
        do.Array(4) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'dualcircular')    
        do.Array(5) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'stacked_circular')    
        do.Array(6) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'spherical')
        do.Array(7) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'circular_rigid_cylinder')    
        do.Array(8) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'stacked_circular_rigid_cylinder')    
        do.Array(9) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'circular_rigid_sphere')    
        do.Array(10) = 1; 
    elseif strcmpi(FigConfig.Arrays{iarray}, 'spherical_rigid_sphere')    
        do.Array(11) = 1; 
    end
end     

array_geo = cell(1,narray_geo);

count=0;
for iarray=1:narray_geo
    array_geo{iarray} = ConfigAllArrays{iarray,1,1}.Array.Geo;
    if do.Array(iarray)
        count=count+1;
        array_geo_disp{count} = ConfigAllArrays{iarray,1,1}.Array.GeoDisp;
    end
end

narrays = sum(do.Array==true);
min_ylim = 1/2/narrays;
nbf=length(methods);
frac = 2/5;

if strcmpi(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'3D') || strcmpi(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'2Daz')
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    if ~isfield(FigConfig.Do,'Subplot') && ~FigConfig.Do.Subplot
                        figure;
                    end
                    ibf = nbf; clear p
                    for imethod=nmethods:-1:1
                        if do.BFMethod(imethod); iarr=1; for iarray=1:narray_geo
                        if do.Array(iarray)
                            
%                             if1 = ifRange(iarray_geo,iM,ir,imethod,isource,inoise,ie,1);
%                             if2 = ifRange(iarray_geo,iM,ir,imethod,isource,inoise,ie,2);
%                             fRange = f(if1:if2);
                            FR = squeeze(BFEvalAllArrays.FRangeaz(iarray,iM,ir,imethod,isource,inoise,ie,:));
%                         semilogx(FR,[min_ylim+(iarray_geo-1)/(narray_geo) min_ylim+(iarray_geo-1)/(narray_geo)], ...
%                             'Color',colours{imethod},'LineStyle',lines{imethod},'LineWidth',1.2); hold on;   
                            
                            if strcmpi(methods{ibf},'LSB')  
                                p(ibf,iarr,1) = semilogx(FR,[min_ylim*(1+frac)+(iarr-1)/(narrays) min_ylim*(1+frac)+(iarr-1)/(narrays)], ...
                                'Color',colours{ibf},'LineStyle',lines{ibf},'LineWidth',lineWidth); hold on; 
                              p(ibf,iarr,2) = semilogx([FRangeLS(iarray,1,1) FRangeLS(iarray,1,1)],[min_ylim*(1)+(iarr-1)/(narrays) min_ylim*(1+2*frac)+(iarr-1)/(narrays)], ...
                                'Color',colours{ibf},'LineStyle','-','LineWidth',lineWidth); hold on;
                             p(ibf,iarr,3) = semilogx([FRangeLS(iarray,1,2) FRangeLS(iarray,1,2)],[min_ylim*(1)+(iarr-1)/(narrays) min_ylim*(1+2*frac)+(iarr-1)/(narrays)], ...
                                'Color',colours{ibf},'LineStyle','-','LineWidth',lineWidth); hold on;
                            elseif strcmpi(methods{ibf},'SDB')  
                                p(ibf,iarr,1) = semilogx(FR,[min_ylim+(iarr-1)/(narrays) min_ylim+(iarr-1)/(narrays)], ...
                                'Color',colours{ibf},'LineStyle',lines{ibf},'LineWidth',lineWidth); hold on; 
                            else
                                p(ibf,iarr,1) = semilogx(FR,[min_ylim*(1-frac)+(iarr-1)/(narrays) min_ylim*(1-frac)+(iarr-1)/(narrays)], ...
                                'Color',colours{ibf},'LineStyle',lines{ibf},'LineWidth',lineWidth); hold on; 
                            end
                            iarr=iarr+1;
                        end
                        
                        end; ibf=ibf-1; end
                    end
                    
                    if FigConfig.Do.Title
                        title('Frequency range in azimuth');
                    end
                    
                    ax = gca;
                    ax.XLabel.Interpreter = 'latex';
                    ax.TickLabelInterpreter = 'latex';
                    if FigConfig.Do.XLabel; xlabel('Frequency (Hz)'); end
                    ax.YLim =[0 1];
                    if FigConfig.Do.YTicks
                        ax.YTick = (0:narrays-1)/narrays + 1/2/narrays;
                        ax.YTickLabel = array_geo_disp;
                    else
                        ax.YTick = '';
                    end
                    ax.YDir = 'reverse';
                    if FigConfig.Do.Legend; legend( p(:,1,1), methods); end
                    if FigConfig.Do.XTicks
                        ax.XTick = foctaves;
                        ax.XTickLabel = foctavesTicks;
                    end
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);
                    grid on;
                    set(gca,'FontSize',FigConfig.FontSize);
                    
                   
                end
            end
        end
    end
end
end



%elevation

if strcmpi(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'3D') || strcmpi(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'2Del')

min_ylim = 1/2/narray_geo;
max_ylim = 1 - 1/2/narray_geo;
imethod=length(methods);
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
            for iM=1:nM
                for ir=1:nr
                    if ~isfield(FigConfig.Do,'Subplot') && ~FigConfig.Do.Subplot
                        figure;
                    end
                    for imethod=nmethods:-1:1
                        if do.BFMethod(imethod); iarr=1; for iarray=1:narray_geo
                        if do.Array(iarray)
                            
%                             if1 = ifRange(iarray_geo,iM,ir,imethod,isource,inoise,ie,1);
%                             if2 = ifRange(iarray_geo,iM,ir,imethod,isource,inoise,ie,2);
%                             fRange = f(if1:if2);
                            FR = squeeze(BFEvalAllArrays.FRangeaz(iarray,iM,ir,imethod,isource,inoise,ie,:));
%                         semilogx(FR,[min_ylim+(iarray_geo-1)/(narray_geo) min_ylim+(iarray_geo-1)/(narray_geo)], ...
%                             'Color',colours{imethod},'LineStyle',lines{imethod},'LineWidth',1.2); hold on;   
                            
                            if strcmpi(methods{imethod},'LSB')  
                                semilogx(FR,[min_ylim*(1+frac)+(iarr-1)/(narrays) min_ylim*(1+frac)+(iarr-1)/(narrays)], ...
                                'Color',colours{imethod},'LineStyle',lines{imethod},'LineWidth',lineWidth); hold on; 
                              semilogx([FRangeLS(iarr,1,1) FRangeLS(iarr,1,1)],[min_ylim*(1)+(iarr-1)/(narrays) min_ylim*(1+2*frac)+(iarr-1)/(narrays)], ...
                                'Color',colours{imethod},'LineStyle','-','LineWidth',lineWidth); hold on;
                              semilogx([FRangeLS(iarr,1,2) FRangeLS(iarr,1,2)],[min_ylim*(1)+(iarr-1)/(narrays) min_ylim*(1+2*frac)+(iarr-1)/(narrays)], ...
                                'Color',colours{imethod},'LineStyle','-','LineWidth',lineWidth); hold on;
                            elseif strcmpi(methods{imethod},'SDB')  
                                semilogx(FR,[min_ylim+(iarr-1)/(narrays) min_ylim+(iarr-1)/(narrays)], ...
                                'Color',colours{imethod},'LineStyle',lines{imethod},'LineWidth',lineWidth); hold on; 
                            else
                                semilogx(FR,[min_ylim*(1-frac)+(iarr-1)/(narrays) min_ylim*(1-frac)+(iarr-1)/(narrays)], ...
                                'Color',colours{imethod},'LineStyle',lines{imethod},'LineWidth',lineWidth); hold on; 
                            end
                            iarr=iarr+1;
                        end
                        
                        end; imethod=imethod-1; end
                    end
                    
                    if FigConfig.Do.Title
                        title('Frequency range in elevation');
                    end
                    
                    ax = gca;
                    ax.XLabel.Interpreter = 'latex';
                    ax.TickLabelInterpreter = 'latex';
                    if FigConfig.Do.XLabel; xlabel('Frequency (Hz)'); end
                    ax.XLim = [500 min(f(end),fbandcutoff(end))];
                    ax.YLim =[0 1];
                     if FigConfig.Do.YTicks
                        ax.YTick = (0:narray_geo-1)/narray_geo + 1/2/narray_geo;
                        ax.YTickLabel = array_geo_disp;
                    else
                        ax.YTick = '';
                     end
                    ax.YDir = 'reverse';
                    if FigConfig.Do.Legend; legend(methods); end
                    if FigConfig.Do.XTicks
                        ax.XTick = foctaves;
                        ax.XTickLabel = foctavesTicks;
                    end
                    xlim([max(f(1),fbandcutoff(1)) min(f(end),fbandcutoff(end))]);
                    grid on;
                    set(gca,'FontSize',FigConfig.FontSize);
                    

                end
            end
        end
    end
end
end

   


end



