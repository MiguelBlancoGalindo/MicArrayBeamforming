function [ha,cb]=PlotDirMapEl(Resp,El,f,FigConfig)
% function [ha,cb]=PlotDirMapEl(Resp,El,f,FigConfig)
% function that plots a given directional response as a function of azimuth
% and frequency. 
%
% input arguments:
%   Resp: directional response as a function of azimuth and frequency. 
%   El: steering elevation vector.
%   f: frequency vector.
%   FigConfig: struct contaning plotting properties.
% output arguments:
%   ha: plot handle.
%   cb: colorbar handle.


load('DirRespColorMap2.mat');
colormap(mycolormap);
fsz=FigConfig.FontSize;
NYticks=8;%not including zero
% [~,El1ind] = min(abs(El(1)-(-89:El(2)-El(1):90)));
% [~,Elendind] = min(abs(El(end)-(-89:El(2)-El(1):90)));

if ~isfield(FigConfig,'NormType')
    FigConfig.NormType='';
end

if strcmp(FigConfig.NormType,'freq')
    %hi=imagesc(El,f,(Resp)-(max(Resp(findnearest(FigConfig.Normf,f),:))),FigConfig.CbarLim); cb=colorbar;
    %ha=imgca;
    ha=surf(El,f,(Resp)-(max(Resp(findnearest(FigConfig.Normf,f),:))),FigConfig.CbarLim); view(2); cb=colorbar;
set(ha,'yscale','log')

elseif strcmp(FigConfig.NormType,'none')
%     hi=imagesc(El,f,Resp,FigConfig.CbarLim); cb=colorbar;
%     ha=imgca;
    ha=surf(El,f,Resp,FigConfig.CbarLim); view(2); cb=colorbar;
    set(ha,'yscale','log')

else
%     hi=imagesc(El,f,(Resp)-max(max(Resp)),FigConfig.CbarLim); cb=colorbar;
%     ha=imgca;
      if isfield(FigConfig,'PlotFunction') && strcmpi(FigConfig.PlotFunction,'pcolor')
          pcolor(Az,f,Resp); 
      else
          surf(Az,f,Resp); view(2)
      end
      %surf(El,f,Resp(:,El1ind:Elendind)); 
      ha=gca;
      
%     h = pcolor(El,f,real(log10(Resp-max(max(Resp)))));
%     h.EdgeColor = 'none';
end
Elticks = -180:90:180;
[~,Elminind] = min(abs(El(1)-Elticks));
[~,Elmaxind] = min(abs(El(end)-Elticks));
Elmin = Elticks(Elminind);
Elmax = Elticks(Elmaxind);
Elticks = Elmin:(Elmax-Elmin)/4:Elmax;

ha.CLim=FigConfig.CbarLim;
if FigConfig.Do.Cbar; cb=colorbar; end
set(gca,'YScale','log')
shading interp
title(FigConfig.Title,'fontsize',fsz)
if FigConfig.Do.XLabel
ha.XLabel.Interpreter = 'latex';
xlabel('$\vartheta \ (^\circ)$','fontsize',fsz);
end
if FigConfig.Do.YLabel
ha.YLabel.Interpreter = 'latex';
ylabel('Frequency (Hz)','fontsize',fsz);
end
ha.XTick = Elticks;
ha.TickLabelInterpreter = 'latex';
set(ha,'YDir','Normal');
xlim([Elmin Elmax]);
ylim([f(1) f(end)]);
foctaves = [63 125 250 500 1000 2000 4000 8000 16000];
foctavesTicks = {'63 ', '125', '250', '500', '1k ', '2k ', '4k ', '8k ', '16k'};
ifirstoctave = find(foctaves>=f(1),1,'first');
ilastoctave = find(foctaves<=f(end),1,'last');
foctaves = foctaves(ifirstoctave:ilastoctave);
foctavesTicks = foctavesTicks(ifirstoctave:ilastoctave);
fcutoff = [63*2^(-1/2),foctaves*2^(1/2)];
ylim([max(f(1),fcutoff(1)) min(f(end),fcutoff(end))]);
if FigConfig.Do.YTicks
set(ha,'YTick',foctaves);
set(ha,'YTickLabel',foctavesTicks);
else
    ha.YTick = '';
end
if exist('cb','var')
set(get(cb,'YLabel'),'Interpreter','latex','string','$20\log_{10}|d(\vartheta,0)|$ \ (dB)','fontsize',fsz);
cb.TickLabelInterpreter = 'latex';
end
if FigConfig.Do.Title;    ha.Title.String = FigConfig.Title; end
set(ha,'FontSize',fsz);
if isfield(FigConfig,'f_max')
hline(FigConfig.f_min/100,'w')
hline(FigConfig.f_max/100,'w')
end

end
