function [ha,cb]=PlotDirMapAz(Resp,Az,f,FigConfig)
% function [ha,cb]=PlotDirMapAz(Resp,Az,f,FigConfig)
% function that plots a given directional response as a function of azimuth
% and frequency. 
%
% input arguments:
%   Resp: directional response as a function of azimuth and frequency. 
%   Az: steering azimuth vector.
%   f: frequency vector.
%   FigConfig: struct contaning plotting properties.
% output arguments:
%   ha: plot handle.
%   cb: colorbar handle.


load('DirRespColorMap2.mat');
colormap(mycolormap);
fsz=FigConfig.FontSize;
NYticks=8;%not including zero

if ~isfield(FigConfig,'NormType')
    FigConfig.NormType='';
end

if strcmp(FigConfig.NormType,'freq')
    %hi=imagesc(Az,f,(Resp)-(max(Resp(findnearest(FigConfig.Normf,f),:))),FigConfig.CbarLim); cb=colorbar;
    %ha=imgca;
    ha=surf(Az,f,(Resp)-(max(Resp(findnearest(FigConfig.Normf,f),:))),FigConfig.CbarLim); view(2); cb=colorbar;
set(ha,'yscale','log')

elseif strcmp(FigConfig.NormType,'none')
%     hi=imagesc(Az,f,Resp,FigConfig.CbarLim); cb=colorbar;
%     ha=imgca;
    ha=surf(Az,f,Resp,FigConfig.CbarLim); view(2); cb=colorbar;
    set(ha,'yscale','log')

else
%     hi=imagesc(Az,f,(Resp)-max(max(Resp)),FigConfig.CbarLim); cb=colorbar;
%     ha=imgca;
      %surf(Az,f,(Resp)-max(max(Resp))); 
      if isfield(FigConfig,'PlotFunction') && strcmpi(FigConfig.PlotFunction,'pcolor')
          pcolor(Az,f,Resp); 
      else
          surf(Az,f,Resp);  view(2)
      end
      ha=gca;
      
%     h = pcolor(Az,f,real(log10(Resp-max(max(Resp)))));
%     h.EdgeColor = 'none';
end

ha.CLim=FigConfig.CbarLim; 
if FigConfig.Do.Cbar; cb=colorbar; end
set(gca,'YScale','log')
shading interp
if FigConfig.Do.XLabel
set(ha.XLabel,'Interpreter','latex','string','$\varphi \ (^\circ)$','fontsize',fsz);
end
ha.XTick=-180:90:180;
ha.TickLabelInterpreter = 'latex';
if FigConfig.Do.YLabel
set(ha.YLabel,'Interpreter','latex','string','Frequency (Hz)','fontsize',fsz);
set(ha,'YDir','Normal');
end
xlim([-180 180]);

foctaves = [63 125 250 500 1000 2000 4000 8000 16000];
foctavesTicks = {'63 ', '125', '250', '500', '1k ', '2k ', '4k ', '8k ', '16k'};
ifirstoctave = find(foctaves>=f(1),1,'first');
ilastoctave = find(foctaves<=f(end),1,'last');
foctaves = foctaves(ifirstoctave:ilastoctave);
foctavesTicks = foctavesTicks(ifirstoctave:ilastoctave);
fcutoff = [63*2^(-1/2),foctaves*2^(1/2)];
ylim([max(f(1),fcutoff(1)) min(f(end),fcutoff(end))]);
if FigConfig.Do.YTicks
    ha.YTick = foctaves;
    ha.YTickLabel=foctavesTicks;
else
    ha.YTick = '';
end
if exist('cb','var') 
    if FigConfig.Do.CbarYLabel
        set(get(cb,'YLabel'),'Interpreter','latex','string','$20\log_{10}|d(0,\varphi)|$ (dB)','fontsize',fsz);
    end
    cb.TickLabelInterpreter = 'latex';
end
if FigConfig.Do.Title;    ha.Title.String = FigConfig.Title; end
set(ha,'FontSize',fsz);
if isfield(FigConfig,'f_max')
hline(FigConfig.f_min/100,'w')
hline(FigConfig.f_max/100,'w')
end

end
