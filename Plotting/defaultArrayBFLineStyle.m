function [FigConfig,ArrayGeoInitials,BFInitials] = defaultArrayBFLineStyle(FigConfig)
% function [FigConfig,ArrayGeoInitials,BFInitials] = defaultArrayBFLineStyle(FigConfig)
% function to obtain the default colours, line styles and acronyms for the
% arrays and beamformers under study. 
%
% input parameters:
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.
%
% output parameters:
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.
%   ArrayGeoInitials: cell array with microphone array geometry acronyms.
%   BFInitials: cell array with beamformer acronyms.


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
    colours = {yellow,orange,red,blue,purple,black,green,blue,black,blue,green};
    lines = {'-','-','-','-','-','-','-','--','--','-.','-.'};
    %for multicircular arrays
    %arrays with same geometry of microphones have same colour (e.g. circular,
    %dual circular, tricircular, rectangular, etc)
    %arrays with the same baffle type have the same line style (e.g. equiangle
    %arrangement, equispacing per circumference or uniform xy grid)
    % colours = {blue,purple,'c',green,purple,'c',green,yellow,black,red,orange,'m'};
    % lines = {'-','-','-','-','--','--','--','-.',':','-.','-.','-.'};
    markers = [];
elseif strcmpi(FigConfig.Compare,'methods')
    colours = {'c','b','r','r',burdeaux,burdeaux,darkGreen,'g'};
    lines = {'-','-','-.','--','-','--','-.','-.'};
    markers = {'none','none','none','none','none','none','none','none'};
    
elseif strcmpi(FigConfig.Compare,'orders')
    colours = [];
    lines = [];
    markers = [];
    
end

FigConfig.LineStyle.Colours = colours;
FigConfig.LineStyle.Lines = lines;
FigConfig.LineStyle.Markers = markers;