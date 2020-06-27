function [derrorFreq,derrorfRange,fRange] = directionalResponseErrorLSB(BFEvalAllArrays, ConfigAllArrays, doPlot)
% function [derrorFreq,derrorfRange,fRange] = directionalResponseErrorLSB(BFEvalAllArrays, ConfigAllArrays, doPlot)
% function that computes the normalised squared error in dB between the
% synthesised responses with the least-squares beamformer (LSB) and the
% equivalent target (desired) responses for each array and directivity
% order.
%
% input arguments: 
%   BFEvalAllArrays: evaluation struct containing several performance
%       metrics for each array, and directivity orders among other things.
%   ConfigAllArrays: cell array containing all the configuration settings 
%       for each microphone array and directivity orders.
%   doPlot: flag to plot the equivalent error for each array and order.
%
% output arguments:
%   derrorFreq: normalised squared error in dB as a function of frequency,
%       for each array and directivity order.
%   derrorfRange: average of derrorFreq across all frequencies within fRange.
%   fRange: frequency range within which derrorFreq is smaller than
%       threshold (specified below as -20 dB).

if isfield(BFEvalAllArrays,'Respaz')
    plane = 'az';
elseif isfield(BFEvalAllArrays,'Respel')
    plane = 'el';
end 
narrays = size(BFEvalAllArrays.(sprintf(strcat('Resp',plane))),1);
nmethods = size(BFEvalAllArrays.(sprintf(strcat('Resp',plane))),4);
nfreqs = size(BFEvalAllArrays.(sprintf(strcat('Resp',plane))),8);
nangles = size(BFEvalAllArrays.(sprintf(strcat('Resp',plane))),9);
nls = 0;
for imethod=1:nmethods; if strcmp(ConfigAllArrays{1,1,1}.Filter.Method{imethod},'ls')
nls = nls+1;
end; end
derrorFreq = zeros(nfreqs,narrays,nls);
derrorfRange = nan(narrays,nls);
f = ConfigAllArrays{1,1,1}.Filter.FreqVector;
fRange=nan(narrays,nls,2);
if1 = cell(narrays,nls);
if2 = cell(narrays,nls);
derror_max = nan(narrays,nls);
derror_min = nan(narrays,nls);
threshold = -20;     %dB
ils = 0;
Config = ConfigAllArrays{1,1,1};
Config.Filter.CurrTarget = Config.Scene.Targets.Angle;
for imethod=1:nmethods; if strcmp(ConfigAllArrays{1,1,1}.Filter.Method{imethod},'ls')
    ils = ils+1;
    
    Config.Filter.Method = ConfigAllArrays{1,1,1}.Filter.Method{imethod};
    if nmethods == length(ConfigAllArrays{1,1,1}.Filter.ls.order)
        Config.Filter.ls.CurrOrder = ConfigAllArrays{1,1,1}.Filter.ls.order(imethod);
        Config.Filter.ls.CurrMode = ConfigAllArrays{1,1,1}.Filter.ls.mode(imethod);
    else
        Config.Filter.ls.CurrOrder = ConfigAllArrays{1,1,1}.Filter.ls.order;
        Config.Filter.ls.CurrMode = ConfigAllArrays{1,1,1}.Filter.ls.mode;
    end
    d_d = getLSIdealPattern(Config).';
for iarray=1:narrays
    %d = 10.^(squeeze(BFEvalAllArrays.(sprintf(strcat('Resp',plane)))(iarray,1,1,imethod,1,1,1,:,:))/20);
    d = squeeze(BFEvalAllArrays.(sprintf(strcat('Resp',plane)))(iarray,1,1,imethod,1,1,1,:,:));
%    derror=abs(d - abs(d_d));
    %sum of "squared pressures"
%    inan = all(isnan(derror),1);
%    derrorFreq(:,iarray,ils)= 10*log10(sum(derror(:,~inan).^2,2)./(sum(abs(d(~inan).').^2)));
    derrorFreq(:,iarray,ils) = directionalResponseError(d,d_d);
    [minerror,iminerror] = min(derrorFreq(:,iarray,ils));
    if1{iarray,ils} = find(derrorFreq(1:iminerror,iarray,ils)<threshold,1,'first');
    if minerror<threshold && ~isempty(if1{iarray,ils})
        if2{iarray,ils}=iminerror + find(derrorFreq(iminerror+1:end,iarray,ils)<threshold,1,'last');
        if minerror<threshold && ~isempty(if2{iarray,ils}) 
            fRange(iarray,ils,1) = f(if1{iarray,ils}); 
            fRange(iarray,ils,2) = f(if2{iarray,ils}); 
        end

    end
    if ~isnan(fRange(iarray,ils,:))
        derrorfRange(iarray,ils) = 10*log10(sum(10.^(derrorFreq(if1{iarray,ils}:if2{iarray,ils},iarray,ils)/10),1)/(if2{iarray,ils}-if1{iarray,ils}+1));
        %derrorfRange(iarray) = 20*log10(sum(10.^(derrorFreq(iffirst:iflast,iarray)/20),1)/(iflast-iffirst+1));
    end
    if ~isempty(if1{iarray,ils}) && ~isempty(if2{iarray,ils}) 
        derror_max(iarray,ils)=max(derrorFreq(if1{iarray,ils}:if2{iarray,ils},iarray,ils));
        derror_min(iarray,ils)=min(derrorFreq(if1{iarray,ils}:if2{iarray,ils},iarray,ils));
    else
        derror_max(iarray,ils)=max(derrorFreq(1:end,iarray,ils));
        derror_min(iarray,ils)=min(derrorFreq(1:end,iarray,ils));
    end
end
[dmin(ils),idmin(ils)] = min(abs(d_d));

end; end

if doPlot
for ils = 1:nls
    FntSz = 13;
    f = ConfigAllArrays{1,1,1}.Filter.FreqVector;
    blue = [0 0.447 0.741];
    orange = [0.85,0.325,0.098];
    red = [1 0 0];
    yellow = [0.929 0.694 0.125];
    purple = [0.494 0.184 0.556];
    green = [0.4660 0.6740 0.188];
    burdeaux = [0.635 0.447 0.741];
    black = [0 0 0];
    ArrayGeoInitials = {'LB','LE','R','C','DC','SC','S','C-RC','SC-RC','C-RS','S-RS'};
    colours = {yellow,orange,red,blue,purple,black,green,blue,black,blue,green};
    lines = {'-','-','-','-','-','-','-','--','--','-.','-.'};

    foctaves = [63 125 250 500 1000 2000 4000 8000 16000];
    foctavesTicks = {'63', '125', '250', '500', '1k', '2k', '4k', '8k', '16k'};
    ifirstoctave = find(foctaves>=f(1),1,'first');
    ilastoctave = find(foctaves<=f(end),1,'last');
    foctaves = foctaves(ifirstoctave:ilastoctave);
    fbandcutoff = [foctaves.*2^(-1/2) foctaves(end).*2^(1/2)];
    
    count=0;
    for iarray=1:narrays
        count=count+1;
        array_geo_disp{count} = ArrayGeoInitials{iarray};
    end
    
    derror_maxi = max(derror_max(:,ils));
    derror_mini = min(derror_min(:,ils));
    max_ylim = derror_maxi + (derror_maxi - derror_mini)/7;
    min_ylim = derror_mini - (derror_maxi - derror_mini)/15;

    figure;
    % limits of f1 and f2
    for iarray=1:narrays
        p1=semilogx([fRange(iarray,ils,1) fRange(iarray,ils,2)], [derror_maxi + (max_ylim-derror_maxi)*iarray/(narrays+1) derror_maxi + (max_ylim-derror_maxi)*iarray/(narrays+1)],...
            'Color',colours{iarray},'LineStyle',lines{iarray}); hold on;
        %,'Parent',ax2
    end
    
    for iarray=1:narrays
        if isnan(fRange(iarray,ils,1)); if1{iarray,ils}=1; end
        if isnan(fRange(iarray,ils,2)); if2{iarray,ils}=length(f);  end
%         p2=semilogx(f(if1{iarray,ils}:if2{iarray,ils}),squeeze(derrorFreq(if1{iarray,ils}:if2{iarray,ils},iarray,ils)),...
%         'Color',colours{iarray},'LineStyle',lines{iarray}); hold on;
        p2=semilogx(f(1:end),squeeze(derrorFreq(1:end,iarray,ils)),...
        'Color',colours{iarray},'LineStyle',lines{iarray}); hold on;
    end
    %legend(array_geo_disp);
    
    ax=gca;
    ax1_pos = ax.Position; % position of first axes
    ax.XTick = foctaves;
    ax.XTickLabel = foctavesTicks;
    ax.FontSize = FntSz;
    ax.XLim = [f(1) min(f(end),fbandcutoff(end))];
    ax.YLim =[min_ylim max_ylim];
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = '<d_{error}(f)> (dB)';
% %     
%     ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
%     ax2.XScale = 'log';
%     ax2.XTick = foctaves;
%     ax2.XTickLabel = foctavesTicks;
%     ax2.YTick = [];
%     ax2.YTickLabel = [];
%     ax2.XLim = [500 min(f(end),fbandcutoff(end))];
%     ax2.YLim =[min_ylim max_ylim];
%     ax2.XLabel.String = 'Frequency (Hz)';
%     ax2.FontSize = FntSz;
    
   
end
end

end



