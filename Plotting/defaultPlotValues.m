function FigConfig=defaultPlotValues(FigConfig)
% function FigConfig=defaultPlotValues(FigConfig)
% function to obtain the default settings to be used when plotting the
% beamformer metrics. 
%
% input and output parameter:
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.

if ~isfield(FigConfig.Do,'NewFig') || isempty(FigConfig.Do.NewFig); FigConfig.Do.NewFig=true; end
if ~isfield(FigConfig.Do,'Subplot') || isempty(FigConfig.Do.Subplot); FigConfig.Do.Subplot=false; end
if ~isfield(FigConfig.Do,'Legend') || isempty(FigConfig.Do.Legend); FigConfig.Do.Legend=true; end
if ~isfield(FigConfig,'FontSize') || isempty(FigConfig.FontSize); FigConfig.FontSize=16; end
if ~isfield(FigConfig,'CbarLim') || isempty(FigConfig.CbarLim); FigConfig.CbarLim=[-30,10]; end
if ~isfield(FigConfig.Do,'CbarYLabel') || isempty(FigConfig.Do.CbarYLabel); FigConfig.Do.CbarYLabel=true; end
if ~isfield(FigConfig,'path') || isempty(FigConfig.Path); FigConfig.Path='~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/'; end
if ~isfield(FigConfig.Do,'XLabel') || isempty(FigConfig.Do.XLabel); FigConfig.Do.XLabel=true; end
if ~isfield(FigConfig.Do,'YLabel') || isempty(FigConfig.Do.YLabel); FigConfig.Do.YLabel=true; end
if ~isfield(FigConfig.Do,'XTicks') || isempty(FigConfig.Do.XTicks); FigConfig.Do.XTicks=true; end
if ~isfield(FigConfig.Do,'YTicks') || isempty(FigConfig.Do.YTicks); FigConfig.Do.YTicks=true; end