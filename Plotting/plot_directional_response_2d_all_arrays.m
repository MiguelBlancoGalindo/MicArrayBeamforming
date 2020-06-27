function plot_directional_response_2d_all_arrays(ConfigAllArrays,ArrayManAllArrays,BFEvalAllArrays, FigConfig)
% function plot_directional_response_2d_all_arrays(ConfigAllArrays,ArrayManAllArrays,BFEvalAllArrays, FigConfig)
% Function that plots the directional response in 2d of the arrays and/or
% beamforming methods under evaluation. 
%
% Note: use PlotDirectionalResponse2DAllArrays.m wherever possible instead, 
% as it is a more up-to-date version. Note the latter uses BFEvalAllArrays 
% as a struct rather than cell arrray thus relying on the function
% rearrange_directional_response_measures.m to rearrange the data.
%
% input arguments:
%   ConfigAllArrays: cell array containing all the configuration settings 
%       for each microphone array, beamforming method, etc.
%   ArrayManAllArrays: cell array containing the array manifold transfer
%       function of the arrays and beamformers under evaluation. 
%   BFEvalAllArrays: cell array containing containing several performance
%       metrics for each array, beamforming method, etc.
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.


FigConfig.FontSize=16;
FigConfig.CbarLim=[-30,10];
narray_geo = size(ConfigAllArrays,1);
nM = size(ConfigAllArrays,2);
nr = size(ConfigAllArrays,3);

array_geo = cell(1,narray_geo);
array_geo_file = cell(1,narray_geo);
M = zeros(1,nM);
r = zeros(1,nr);

%temporary fix until all config files being used have the Geo
if isfield(ConfigAllArrays{1,1,1},'Array')
    ConfigArray = ConfigAllArrays{1,1,1}.Array;
    if ~isfield(ConfigArray,'Geo')
        arraygeoexists=0;
        array_geo_file={'linear','rectangular','circular','dualcircular','spherical',...
        'circular_rigid_cylinder','circular_rigid_sphere','spherical_rigid_sphere'};
    else
        arraygeoexists=1;
    end
end
for iarray_geo=1:narray_geo
    array_geo{iarray_geo} = ConfigAllArrays{iarray_geo,1,1}.Array.GeoDisp;
    if arraygeoexists==1
        array_geo_file{iarray_geo} = ConfigAllArrays{iarray_geo,1,1}.Array.Geo;
    end
end

for iM=1:nM
    M(iM) = size(ConfigAllArrays{1,iM,1}.Array.MicPos,1);
end
%radius estimated from half the aperture of linear array
for ir=1:nr
    r(ir)= round(ConfigAllArrays{1,1,ir}.Array.Aperture/2,2);
end  
rmm = r*1000;   %radius in mm

nsources=size(BFEvalAllArrays{1,1,1},2);
nnoises=size(BFEvalAllArrays{1,1,1}{1,1},2);
nmethods = size(ConfigAllArrays{1,1,1}.Filter.Method,2);
source1 = 1;
sourcen = 1;
noise1 = 1;
noisen = nnoises;

nmetrics = 6; 
ntotalmethods = 6;
do.metric = zeros(nmetrics);
do.bfmethod = zeros(ntotalmethods);
do.array = zeros(narray_geo);

if ~isfield(FigConfig,'NoiseType') || isempty(FigConfig.NoiseType)
    FigConfig.NoiseType=cell(1,nnoises);
else
    NoiseType=FigConfig.NoiseType(1:end);
    FigConfig.NoiseType=cell(1,nnoises);
    FigConfig.NoiseType(2:end)=strcat('_',NoiseType);
    FigConfig.NoiseType{1} = '_ideal';
end

for imethod=1:length(FigConfig.methods)
    if strcmp(FigConfig.methods{imethod}, 'all')
        do.bfmethod(:) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'ds')
        do.bfmethod(1) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'sda')
        do.bfmethod(2) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'mvdr')
        do.bfmethod(3) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'lcmv')    
        do.bfmethod(4) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'ls')    
        do.bfmethod(5) = 1;
    elseif strcmp(FigConfig.methods{imethod}, 'acc')    
        do.bfmethod(6) = 1;
    end
end 

for iarray_geo=1:length(FigConfig.arrays)
    if strcmp(FigConfig.arrays{iarray_geo}, 'all')
        do.array(:) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'linear_bs')
        do.array(1) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'linear_ef')
        do.array(2) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'rectangular')
        do.array(3) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'circular')    
        do.array(4) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'dualcircular')    
        do.array(5) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'stacked_circular')    
        do.array(6) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'spherical')
        do.array(7) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'circular_rigid_cylinder')    
        do.array(8) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'stacked_circular_rigid_cylinder')    
        do.array(9) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'circular_rigid_sphere')    
        do.array(10) = 1; 
    elseif strcmp(FigConfig.arrays{iarray_geo}, 'spherical_rigid_sphere')    
        do.array(11) = 1; 
    end
end 

if strcmp(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'3D') || strcmp(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'2Daz')

for inoise=noise1:noisen
    for imethod=1:nmethods
        if do.bfmethod(imethod); for isource=source1:sourcen
            for iM=1:nM
                for ir=1:nr
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        figure;
                        FigConfig.Title=[ConfigAllArrays{iarray_geo,iM,ir}.Filter.Method{1,imethod} ' ' array_geo{iarray_geo} ' M=' num2str(M(iM)) ' r=' num2str(r(ir))];
                        if FigConfig.Save; FigConfig.Title=''; end
                        PlotDirMapAz(BFEvalAllArrays{iarray_geo,iM,ir}{imethod,isource}(inoise).Respaz,ArrayManAllArrays{iarray_geo,iM,ir}.Monitor(1).LookAng.Az,ArrayManAllArrays{iarray_geo,iM,ir}.Monitor(1).f,FigConfig);
                        %PlotDirMap_MB(BFEvalAllArrays{o,m}(nnoise).Resp,ArrayManAllArrays{o,m}.Monitor(1).LookAng,ArrayManAllArrays{o,m}.Monitor(1).f,FigConfig);
                        %save(gcf,['Directional_response_' Config.Filter.Method{1,1} '_' array_geo{o} '_M' num2str(M) '_' Config.Array.LookDistance '.pdf']);
                        %save figure
                        hold on;
                        plot([rad2deg(ArrayManAllArrays{iarray_geo,iM,ir}.Monitor(1).LookAng.Az(1)) ...
                            rad2deg(ArrayManAllArrays{iarray_geo,iM,ir}.Monitor(1).LookAng.Az(end))],[ConfigAllArrays{iarray_geo,iM,ir}.Array.falias ConfigAllArrays{iarray_geo,iM,ir}.Array.falias],'--r','LineWidth',2);
                        method = ConfigAllArrays{iarray_geo,iM,ir}.Filter.Method{1,imethod};

                        if strcmp(method,'ls')
                            order = ['_n' num2str(ConfigAllArrays{iarray_geo,iM,ir}.Filter.ls.order) '_'];
                            pattern = ConfigAllArrays{iarray_geo,iM,ir}.Filter.ls.mode;
                            if strcmp(pattern,'hypercardioid')
                                pattern = 'hyp';
                            end
                        else
                            order ='';
                            pattern = '';
                        end
%                         if nnoises>1 && inoise==1
%                             SaveCurr=0;
%                         else
%                             SaveCurr=FigConfig.Save;
%                         end
%                         if SaveCurr==1
                        if FigConfig.Save==1
                            print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                                'Directional_response_fullbw_az_' method order pattern '_' array_geo_file{iarray_geo} '_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{iarray_geo,iM,ir}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                        end
                        end
                    end

                end
            end
        end; end
    end
end
end

% Elevation
if strcmp(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'3D') || strcmp(ConfigAllArrays{1,1,1}.ArrayMan.SteerSpace,'2Del')

for inoise=noise1:noisen
    for imethod=1:nmethods
        if do.bfmethod(imethod); for isource=source1:sourcen
            for iM=1:nM
                for ir=1:nr
                    for iarray_geo=1:narray_geo
                        if do.array(iarray_geo)
                        Respel = BFEvalAllArrays{iarray_geo,iM,ir}{imethod,isource}(inoise).Respel;
                        Angles = ArrayManAllArrays{iarray_geo,iM,ir}.Monitor(1).LookAng.El;
                        figure;
                        FigConfig.Title=[ConfigAllArrays{iarray_geo,iM,ir}.Filter.Method{1,imethod} ' ' array_geo{iarray_geo} ' M=' num2str(M(iM)) ' r=' num2str(r(ir))];              
                        if FigConfig.Save; FigConfig.Title=''; end
                        PlotDirMapEl(Respel,Angles,ArrayManAllArrays{iarray_geo,iM,ir}.Monitor(1).f,FigConfig);
                        %PlotDirMap_MB(BFEvalAllArrays{o,m}(nnoise).Resp,ArrayManAllArrays{o,m}.Monitor(1).LookAng,ArrayManAllArrays{o,m}.Monitor(1).f,FigConfig);
                        %save(gcf,['Directional_response_' Config.Filter.Method{1,1} '_' array_geo{o} '_M' num2str(M) '_' Config.Array.LookDistance '.pdf']);
                        %save figure
                        hold on;
                        plot([rad2deg(ArrayManAllArrays{iarray_geo,iM,ir}.Monitor(1).LookAng.El(1)) ...
                            rad2deg(ArrayManAllArrays{iarray_geo,iM,ir}.Monitor(1).LookAng.El(end))],[ConfigAllArrays{iarray_geo,iM,ir}.Array.falias ConfigAllArrays{iarray_geo,iM,ir}.Array.falias],'--r','LineWidth',2);

                        method = ConfigAllArrays{iarray_geo,iM,ir}.Filter.Method{1,imethod};
                        if strcmp(method,'ls')
                                order = ['_n' num2str(ConfigAllArrays{iarray_geo,iM,ir}.Filter.ls.order) '_'];
                                pattern = ConfigAllArrays{iarray_geo,iM,ir}.Filter.ls.mode;
                                if strcmp(pattern,'hypercardioid')
                                    pattern = 'hyp';
                                end
                            else
                            order ='';
                            pattern = '';
                        end
%                         if nnoises>1 && inoise==1
%                             SaveCurr=0;
%                         else
%                             SaveCurr=FigConfig.Save;
%                         end
%                         if SaveCurr==1
                        if FigConfig.Save==1
                            print(gcf, '-depsc', ['~/Dropbox/Surrey Uni/PhD/Matlab/Beamformers/Phils_toolbox/images/' ...
                                'Directional_response_fullbw_el_' method order pattern '_' array_geo_file{iarray_geo} '_M' num2str(M(iM)) '_r' num2str(rmm(ir)) 'mm_' ConfigAllArrays{iarray_geo,iM,ir}.ArrayMan.LookDistance FigConfig.NoiseType{inoise} '.eps']);
                        end
                        end
                    end

                end
            end
        end; end
    end
end


end