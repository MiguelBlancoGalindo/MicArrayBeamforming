function PlotDirectionalResponse3DAllArrays(ConfigAllArrays,BFEvalAllArrays, FigConfig)
% function PlotDirectionalResponse3DAllArrays(ConfigAllArrays,BFEvalAllArrays, FigConfig)
% Function that plots the directional response in 3d of the arrays and/or
% beamforming methods under evaluation. 
%
% input arguments:
%   ConfigAllArrays: cell array containing all the configuration settings 
%       for each microphone array, beamforming method, etc.
%   BFEvalAllArrays: struct containing several performance metrics for each 
%       array, beamforming method, etc. Note BFEvalAllArrays is a struct 
%       from the output of rearrange_directional_response_measures.m. 
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.


if ~strcmpi(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'3D') 
    error('3D plots not valid for 2D SteerSpaces');
else
    
if ~isfield(FigConfig,'FontSize') || isempty(FigConfig.FontSize); FigConfig.FontSize=8; end
FigConfig=defaultPlotValues(FigConfig);

f=ConfigAllArrays{1,1,1}.Filter.FreqVector;
foctaves = [63 125 250 500 1000 2000 4000 8000 16000];
foctavesTicks = {'63 ', '125', '250', '500', '1k ', '2k ', '4k ', '8k ', '16k'};

phi = ConfigAllArrays{1,1,1}.ArrayMan.LookAng.Az;
theta = ConfigAllArrays{1,1,1}.ArrayMan.LookAng.El;
theta_t = ConfigAllArrays{1,1,1}.Scene.Targets.Angle(:,1);
phi_t = ConfigAllArrays{1,1,1}.Scene.Targets.Angle(:,2);
nTars = length(phi_t);
itar = zeros(nTars,1);
for it = 1:nTars
    [~,~,itar(it)] = getNearestAngle(ConfigAllArrays{1},theta_t(it),phi_t(it));
end

plane = '3d';
BFInitials = {'DSB','SDB','MVDRB1','MVDRB2','LCMVB1','LCMVB2','LSB'};

%method to include part of a string as a field of a struct by using sprintf
%of the concatenated string ("DI"=part of variable + "plane"=newstring) and
%its output (string) include it in brackets to convert it to a field
narrays = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),1);
nM = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),2);
nr = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),3);
nmethods = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),4);
nsources = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),5);
nnoises = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),6);
nepsilons = size(BFEvalAllArrays.(sprintf(strcat('DI',plane))),7);
nfreqs = length(f);

source1 = 1;
sourcen = 1;
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
rmm = r*1000;   %radius in mm

ntotalmethods = nmethods;
do.bfmethod = zeros(ntotalmethods,1);
do.array = zeros(narrays,1);
do.noise = zeros(nnoises,1);
for inoise=1:length(FigConfig.Noises)
    if FigConfig.Noises(inoise)<=nnoises
        do.noise(FigConfig.Noises(inoise)) = true;
    else
        error('specified measurement noise index to be plotted exceeds maximum index');
    end
end

if ~isfield(FigConfig,'NoiseType') || isempty(FigConfig.NoiseType)
    FigConfig.NoiseType=cell(1,nnoises);
else
    noisen = length(FigConfig.Noises);
    if noisen~=length(FigConfig.NoiseType)
        error(['Number of ''' ,'NoiseType', ''' must match number of ''','Noises','''']);
%     else
%         if noisen>1
%             NoiseType=cell(1,noisen);
%             NoiseType(2:end)=FigConfig.NoiseType(1:end);
%             FigConfig.NoiseType=strcat('_',NoiseType);
%         elseif noisen==1; FigConfig.NoiseType=strcat('_',FigConfig.NoiseType);
%         end
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

for iarr=1:length(FigConfig.Arrays)
    if strcmpi(FigConfig.Arrays{1}, 'all')
        do.array(:) = 1;
        break
    else
        for iarray = 1:narrays
            if strcmpi(ConfigAllArrays{iarray,1,1}.Array.Geo, FigConfig.Arrays{iarr})
                do.array(iarray) = 1;
            end
        end
    end
end 
array_geo = cell(1,narrays);
array_geo_disp = cell(1,narrays);

for iarray=1:narrays
    array_geo{iarray} = ConfigAllArrays{iarray,1,1}.Array.Geo;
    array_geo_disp{iarray} = ConfigAllArrays{iarray,1,1}.Array.GeoDisp;
end

if isfield(FigConfig,'Freq') && ~isempty(FigConfig.Freq) && ...
    (~isfield(FigConfig,'OctBand') || ~isempty(FigConfig.OctBand))
    avg = 0;
    fplot = FigConfig.Freq;
    nfreqs = length(FigConfig.Freq);
    
elseif isfield(FigConfig,'OctBand') && ~isempty(FigConfig.OctBand) && ...
    (~isfield(FigConfig,'Freq') || ~isempty(FigConfig.Freq))
    avg=1;
    fplot = FigConfig.OctBand;
    nfreqs = length(FigConfig.OctBand);
    fdiff = zeros(nfreqs,1);
    for ifreq=1:nfreqs
        fdiff(ifreq) = min(foctaves-fplot(ifreq));
    end
    if sum(fdiff)>0
        error('requested frequencies do not match standard octave band centre frequencies');
    end
else
    error(['Select only values for ''', 'Freq', ''' or ''', 'OctBand', '''']);
end

ifreq = zeros(nfreqs,1);
icutoff = zeros(nfreqs,2);
fplotprint = cell(nfreqs,1);
for kf=1:nfreqs
    [~,ifreq(kf)] = min(abs(fplot(kf)-f));
    fcutoff = [fplot(kf)*2^(-1/2),fplot(kf)*2^(1/2)];
    [~,icutoff(kf,1)] = min(abs(fcutoff(1)-f));
    [~,icutoff(kf,2)] = min(abs(fcutoff(2)-f)); 
    icutoff(kf,2) = icutoff(kf,2) -1;
    if fplot(kf)>=1000
        fplotprint{kf} = strcat(num2str(fplot(kf)/1000),' kHz');
    else
        fplotprint{kf} = strcat(num2str(fplot(kf)),' Hz');
    end
end
nplots = nfreqs;
doMag = false; doPhase = false; doAz = false; doEl = false;
if isfield(FigConfig.Do.OnAxis,'Mag') && FigConfig.Do.OnAxis.Mag
    nplots = nplots+2;
    doMag = true;
end
if isfield(FigConfig.Do.OnAxis,'Phase') && FigConfig.Do.OnAxis.Phase
    nplots = nplots+2;
    doPhase = true;
end
if isfield(FigConfig.Do,'Az') && FigConfig.Do.Az
    doAz = true;
end
if isfield(FigConfig.Do,'El') && FigConfig.Do.El
    doEl = true;
end

FigConfig.Units = 'points';
colwidth = 252; % pt
heightscale = 0.38*ceil(nplots/2);

% for 3D plots
FigConfig.Subplot.Width =  0.38 * colwidth;
FigConfig.Subplot.Height = 1/ceil(nplots/2)*0.9 * heightscale * colwidth;   
FigConfig.Subplot.Gap = 0.08 * colwidth;
FigConfig.Subplot.Pos.Left = 0.06* colwidth + FigConfig.Subplot.Width*(0:1) + FigConfig.Subplot.Gap*2*(0:1);
FigConfig.Subplot.Pos.Bottom = 15 + FigConfig.Subplot.Height*(ceil(nplots/2)-1:-1:0) + FigConfig.Subplot.Gap*(ceil(nplots/2)-1:-1:0);

%for on axis plots
FigConfig.Subplot2.Width =  0.33 * colwidth;
FigConfig.Subplot2.Height = 1/ceil(nplots/2)*0.5 * heightscale * colwidth;   
FigConfig.Subplot2.Gap = 0.05 * colwidth;
FigConfig.Subplot2.Pos.Left = 0.1* colwidth + FigConfig.Subplot2.Width*(0:1) + FigConfig.Subplot2.Gap*2*(0:1);
FigConfig.Subplot2.Pos.Bottom = 15 + [(FigConfig.Subplot.Height+FigConfig.Subplot.Gap)*(ceil(nfreqs/2)) + (FigConfig.Subplot2.Height+FigConfig.Subplot2.Gap)*(ceil((nplots-nfreqs)/2)-1:-1:0) FigConfig.Subplot.Height*(ceil(nfreqs/2)-1:-1:0)] +15;

%for 2D plot
FigConfig.Subplot1.Width =  0.48 * colwidth;
FigConfig.Subplot1.Height = FigConfig.Subplot1.Width *0.6;   
FigConfig.Subplot1.Gap = 0.01 * colwidth;
FigConfig.Subplot1.Pos.Left = 0.03* colwidth + FigConfig.Subplot1.Width*(0:1) + FigConfig.Subplot1.Gap*(0:1);
FigConfig.Subplot1.Pos.Bottom = FigConfig.Subplot2.Pos.Bottom(1)+FigConfig.Subplot2.Height-FigConfig.Subplot1.Height;

for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=1:nsources; in=0;
        for inoise=1:nnoises; if do.noise(inoise); in=in+1;
            for iM=1:nM
                for ir=1:nr
                    for iarray=1:narrays
                        if do.array(iarray)
                            fig=figure('Units',FigConfig.Units,'Position',[0,0,colwidth,heightscale*colwidth]);
                            Resp = squeeze(BFEvalAllArrays.Resp3d(iarray,iM,ir,imethod,isource,inoise,ie,:,:));
                            Resp_t = squeeze(Resp(:,itar(isource)));
                            %removing the linear phase from the response
                            phase = unwrap(angle(Resp_t));  
                            [~,if1] = min(abs(f-20));             
                            [~,if2] = min(abs(f-20000));         
                            %calculate group delay to remove linear phase
                            d = f(if1:if2)'\phase(if1:if2);         
                            phase_min = phase - d*f.';      
                            %normalise the response at 1kHz
                            [~,if1k] = min(abs(f-1000));
                            phase_1k = phase_min(if1k);
                            phase_min = phase_min - phase_1k;  
                            %Wrap between -pi e +pi
                            phase_min = 180/pi*atan2(sin(phase_min),cos(phase_min));    
                            iplot=0;
                            if doMag
                                iplot=iplot+1;
                                ax=axes('Units',FigConfig.Units,'Position',[FigConfig.Subplot2.Pos.Left(1),FigConfig.Subplot2.Pos.Bottom(1),FigConfig.Subplot2.Width,FigConfig.Subplot2.Height]);
                                semilogx(f,db(abs(Resp_t)));
                                ax.XLabel.Interpreter = 'latex';
                                ax.YLabel.Interpreter = 'latex';
                                ax.TickLabelInterpreter = 'latex';
                                xlabel('f (Hz)');
                                ylabel('$|d(\Omega_l)|^2$ (dB)'); 
                                xlim([f(1) f(end)]);
                                ylim([-20 +10]);
                                set(gca,'XTick',[100,1000,10000]);
                                set(gca,'XTickLAbel',{'100','1k','10k'});
                                set(gca,'YTick',[-20,-10,0,10]);
                                set(gca,'YTickLabel',{'-20','-10','0','10'});
                                set(gca,'XGrid','on');
                                set(gca,'YGrid','on');
                                title(array_geo_disp(iarray),'Interpreter','latex','fontsize',FigConfig.FontSize);
                                titleDone = 1;
                            end
                            if doPhase
                                iplot=iplot+1;
                                ax=axes('Units',FigConfig.Units,'Position',[FigConfig.Subplot2.Pos.Left(1),FigConfig.Subplot2.Pos.Bottom(2),FigConfig.Subplot2.Width,FigConfig.Subplot2.Height]);
                                semilogx(f,phase_min)
                                ax.XLabel.Interpreter = 'latex';
                                ax.YLabel.Interpreter = 'latex';
                                ax.TickLabelInterpreter = 'latex';
                                xlabel('f (Hz)');
                                ylabel('$\angle d(\Omega_l)(^\circ)$');
                                xlim([f(1) f(end)]);
                                ylim([-180 180]);
                                set(gca,'XTick',[100,1000,10000]);
                                set(gca,'XTickLabel',{'100','1k','10k'});
                                set(gca,'YTick',[-180,-90,0,90,180]);
                                ax.YLabel.Position(1) = 8.37;
                                %set(gca,'YTickLabel',{'-$\pi$','-$\pi$/2','0','$\pi$/2','$\pi$'});
                                set(gca,'XGrid','on');
                                set(gca,'YGrid','on');
                                if ~titleDone
                                    title(array_geo_disp(iarray),'Interpreter','latex','fontsize',FigConfig.FontSize);
                                    titleDone = 1;
                                end
                                    
                            end

                            if doAz && doEl
                                Az = phi(theta==0);
                                iAz = find(theta==0);
                                [Az,ind] = sort(Az);
                                Respaz = db(Resp(:,iAz(ind)));
                                ax=axes('Units',FigConfig.Units,'Position',[FigConfig.Subplot2.Pos.Left(2),FigConfig.Subplot2.Pos.Bottom(1),FigConfig.Subplot2.Width,FigConfig.Subplot2.Height]);
                                FigConfig.Do.Cbar=true;
                                PlotDirMapAz(Respaz,Az,f,FigConfig);

                                El = theta(phi==0);
                                iEl = find(phi==0);
                                [El,ind] = sort(El);
                                Respel = db(Resp(:,iEl(ind)));
                                ax=axes('Units',FigConfig.Units,'Position',[FigConfig.Subplot2.Pos.Left(2),FigConfig.Subplot2.Pos.Bottom(2),FigConfig.Subplot2.Width,FigConfig.Subplot2.Height]);
                                PlotDirMapEl(Respel,El,f,FigConfig);
                           
                            else
                                if doAz
                                Az = phi(theta==0);
                                iAz = find(theta==0);
                                [Az,ind] = sort(Az);
                                Respaz = db(Resp(:,iAz(ind)));
                                ax=axes('Units',FigConfig.Units,'Position',[FigConfig.Subplot1.Pos.Left(2),FigConfig.Subplot1.Pos.Bottom,FigConfig.Subplot1.Width,FigConfig.Subplot1.Height]);
                                FigConfig.Do.Cbar=true;
                                [~,cb]=PlotDirMapAz(Respaz,Az,f,FigConfig);
                                set(get(cb,'YLabel'),'Interpreter','latex','string','$|d(\Omega)|^2$ (dB)');
                                
                                elseif doEl
                                El = theta(phi==0);
                                iEl = find(phi==0);
                                [El,ind] = sort(El);
                                Respel = db(Resp(:,iEl(ind)));
                                ax=axes('Units',FigConfig.Units,'Position',[FigConfig.Subplot1.Pos.Left(2),FigConfig.Subplot1.Pos.Bottom,FigConfig.Subplot1.Width,FigConfig.Subplot1.Height]);
                                PlotDirMapEl(Respel,El,f,FigConfig);
                                end
                            end

                            iplot = iplot*2;

                            for kf=1:nfreqs
                                if avg
                                    %phase response at centre frequency
                                    phase = angle(Resp(icutoff(kf,1),:) - d*fplot(kf) - phase_1k);
                                    %mg response as average across band
                                    Respf = mean( abs(Resp(icutoff(kf,1):icutoff(kf,2),:)) .* exp(1i*phase));
                                else
                                    Respf = Resp(ifreq(kf),:);    
                                end
                                iplot = iplot+1;
                                icol = mod(iplot,2); if icol==0; icol=2; end
                                irow = ceil(iplot/2);
                                ax=axes('Units',FigConfig.Units,'Position',[FigConfig.Subplot.Pos.Left(icol),FigConfig.Subplot.Pos.Bottom(irow),FigConfig.Subplot.Width,FigConfig.Subplot.Height]);
                                FigConfig.Title=[ConfigAllArrays{iarray,iM,ir}.Filter.Method{imethod} ' ' array_geo_disp{iarray} ' ' FigConfig.NoiseType{in}];
                                if FigConfig.Do.Save; FigConfig.Title=''; FigConfig.Do.Cbar=false; else FigConfig.Do.Cbar=true; end
                                polar3d(Respf,[theta phi],FigConfig);
                                title(sprintf([fplotprint{kf} ' (%0.1f dB)'],max(db(Respf))),'Interpreter','latex','fontsize',FigConfig.FontSize);
                                ax.Title.Position(3) = ax.Title.Position(3)*0.4;
                                method = ConfigAllArrays{iarray,iM,ir}.Filter.Method{imethod};
                                if strcmpi(method,'ls')
                                    order = ['_n' num2str(ConfigAllArrays{iarray,iM,ir}.Filter.ls.order) '_'];
                                    pattern = ConfigAllArrays{iarray,iM,ir}.Filter.ls.mode;
                                    if strcmpi(pattern,'hypercardioid')
                                        pattern = 'hyp';
                                    end
                                else
                                    order ='';
                                    pattern = '';
                                end
                            end
                            cb1 = colorbar;
                            cb1.TickLabelInterpreter = 'latex';
                            cb1.Ticks = [0 45 90 135 180];
                            set(get(cb1,'YLabel'),'Interpreter','latex','string','$\angle d(\Omega)(^\circ)$');
                            cb1.Position = [0.5315    0.0303    0.0270    0.1398];

                            if FigConfig.Do.Save==1
                                print(gcf, [FigConfig.Path ...
                                    'Directional_response_fullbw_az_' method order pattern '_' array_geo{iarray} '_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{iarray,iM,ir}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps'], '-depsc','-opengl');
    %                             print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
    %                                 'Directional_response_fullbw_az_' method order pattern '_' array_geo{iarray_geo} '_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{iarray_geo,iM,ir}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                            end
                            
                        end
                    end

                end
            end
        end;  end
    end; end
end
end



end


end