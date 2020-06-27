function PlotDirectionalResponse2DAllArrays(ConfigAllArrays,BFEvalAllArrays, FigConfig)
% function PlotDirectionalResponse2DAllArrays(ConfigAllArrays,BFEvalAllArrays, FigConfig)
% Function that plots the directional response in 2d of the arrays and/or
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


FigConfig=defaultPlotValues(FigConfig);

f=ConfigAllArrays{1,1,1}.Filter.FreqVector;
%Rather than loading ArrayMan which is heavy, use Config to get angle
%vector

Az = sort(ConfigAllArrays{1,1,1}.ArrayMan.LookAng.Az(ConfigAllArrays{1,1,1}.ArrayMan.LookAng.El==0));
El1 = ConfigAllArrays{1,1,1}.ArrayMan.LookAng.El(ConfigAllArrays{1,1,1}.ArrayMan.LookAng.Az==0);
El2 = -ConfigAllArrays{1,1,1}.ArrayMan.LookAng.El(ConfigAllArrays{1,1,1}.ArrayMan.LookAng.Az==180) + 180;
El2(El2>180) = El2(El2>180)-360;
El = sort([El1; El2]);

if isfield(BFEvalAllArrays,'DIaz')
    plane = 'az';
elseif isfield(BFEvalAllArrays,'DIel')
    plane = 'el';
end

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

if strcmpi(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'3D') || strcmpi(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'2Daz')

    
for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen; in=0;
        for inoise=1:nnoises; if do.noise(inoise); in=in+1;
            for iM=1:nM
                for ir=1:nr
                    for iarray=1:narrays
                        if do.array(iarray)
                        Respaz = squeeze(BFEvalAllArrays.Respaz(iarray,iM,ir,imethod,isource,inoise,ie,:,:));
                        figure;
                        FigConfig.Title=[ConfigAllArrays{iarray,iM,ir}.Filter.Method{imethod} ' ' array_geo_disp{iarray} ' ' FigConfig.NoiseType{in}];
                        if FigConfig.Do.Save; FigConfig.Title=''; FigConfig.Do.Cbar=false; else FigConfig.Do.Cbar=true; end
                        PlotDirMapAz(db(Respaz),Az,f,FigConfig);
                        hold on;
                        plot([Az(1) ...
                            Az(end)],[ConfigAllArrays{iarray,iM,ir}.Array.falias ConfigAllArrays{iarray,iM,ir}.Array.falias],'--r','LineWidth',2);
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

                        if FigConfig.Do.Save==1
                            print(gcf, [FigConfig.Path ...
                                'Directional_response_fullbw_az_' method order pattern '_' array_geo{iarray} '_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{iarray,iM,ir}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps'], '-depsc','-opengl');
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

% Elevation
if strcmpi(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'3D') || strcmpi(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'2Del')

for ie=epsilon1:epsilonn
for imethod=1:nmethods
    if do.bfmethod(imethod); for isource=source1:sourcen
        for inoise=1:nnoises; if do.noise(inoise)
            for iM=1:nM
                for ir=1:nr
                    for iarray=1:narrays
                        if do.array(iarray)
                        Respel = squeeze(BFEvalAllArrays.Respel(iarray,iM,ir,imethod,isource,inoise,ie,:,:));
                        figure;
                        FigConfig.Title=[ConfigAllArrays{iarray,iM,ir}.Filter.Method{imethod} ' ' array_geo_disp{iarray} ' ' FigConfig.NoiseType{inoise}];              
                        if FigConfig.Do.Save; FigConfig.Title=''; FigConfig.Do.Cbar=false; else FigConfig.Do.Cbar=true; end
                        PlotDirMapEl(db(Respel),El,f,FigConfig);
                        hold on;
                        plot([El(1) ...
                            El(end)],[ConfigAllArrays{iarray,iM,ir}.Array.falias ConfigAllArrays{iarray,iM,ir}.Array.falias],'--r','LineWidth',2);

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

                        if FigConfig.Do.Save==1                          
                            print(gcf, [FigConfig.Path ...
                                'Directional_response_fullbw_el_' method order pattern '_' array_geo{iarray} '_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{iarray,iM,ir}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps'], '-depsc','-opengl');

                        end
                        end
                    end

                end
            end
        end; end
    end; end
end
end

end

end