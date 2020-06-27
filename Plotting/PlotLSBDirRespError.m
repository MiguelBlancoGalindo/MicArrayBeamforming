function PlotLSBDirRespError(BFEvalAllArrays,ConfigAllArrays, derror, FigConfig)
% function PlotLSBDirRespError(BFEvalAllArrays,ConfigAllArrays, derror, FigConfig)
% Function that plots the error in the synthesised directivity response of 
% the least-squares beamformer (LSB) compared to the equivalent target 
% response.
%
% input arguments:
%   BFEvalAllArrays: struct containing several performance metrics for each 
%       array, beamforming method, etc. Note BFEvalAllArrays is a struct 
%       from the output of rearrange_directional_response_measures.m. 
%   ConfigAllArrays: cell array containing all the configuration settings 
%       for each microphone array, beamforming method, etc.
%   derror: normalised squared error in dB between syntehsised and target
%       directional responses.
%   FigConfig: struct contaning properties of which metrics, arrays and
%       beamformers to plot, among all included in BFEvalAllArrays.


if ~isfield(FigConfig,'FontSize') || isempty(FigConfig.FontSize); FigConfig.FontSize = 13; end
f = ConfigAllArrays{1,1,1}.Filter.FreqVector;

ArrayGeoInitials = {'LB','LE','R','C','DC','SC','S','C-RC','SC-RC','C-RS','S-RS'};

if isfield(BFEvalAllArrays,'DIaz')
    plane = 'az';
elseif isfield(BFEvalAllArrays,'DIel')
    plane = 'el';
end

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
foctaves = [63 125 250 500 1000 2000 4000 8000 16000];
foctavesTicks = {'63', '125', '250', '500', '1k', '2k', '4k', '8k', '16k'};
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

if isfield(FigConfig,'Type') && strcmp(FigConfig.Type,'difference') 
    FigConfig.Error = 'error';
else
    FigConfig.Error = '';
end
    
rmm=r*1000;
nmetrics = 9; 
ntotalmethods = nmethods;
do.array = zeros(narrays,1);


for iarray_geo=1:length(FigConfig.Arrays)
    if strcmp(FigConfig.Arrays{iarray_geo}, 'all')
        do.array(:) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'linear_bs')
        do.array(1) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'linear_ef')
        do.array(2) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'rectangular')
        do.array(3) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'circular')    
        do.array(4) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'dualcircular')    
        do.array(5) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'stacked_circular')    
        do.array(6) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'spherical')
        do.array(7) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'circular_rigid_cylinder')    
        do.array(8) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'stacked_circular_rigid_cylinder')    
        do.array(9) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'circular_rigid_sphere')    
        do.array(10) = 1; 
    elseif strcmp(FigConfig.Arrays{iarray_geo}, 'spherical_rigid_sphere')    
        do.array(11) = 1; 
    end
end     

array_geo = cell(1,narrays);

count=0;
for iarray_geo=1:narrays
    array_geo{iarray_geo} = ConfigAllArrays{iarray_geo,1,1}.Array.Geo;
    if do.array(iarray_geo)
        count=count+1;
        %array_geo_disp{count} = ConfigAllArrays{iarray_geo,1,1}.Array.GeoDisp;
        array_geo_disp{count} = ArrayGeoInitials{iarray_geo};
    end
end

narrays = length(array_geo_disp);
min_ylim = 1/2/narrays;
max_ylim = 1 - 1/2/narrays;
derror1 = zeros(length(f),narrays+1);
derror1(:,1:narrays) = derror(:,do.array==1);
for ie=epsilon1:epsilonn
    for isource=source1:sourcen
        for inoise=noise1:noisen
           
                    if ~isfield(FigConfig.Do,'Subplot') && ~FigConfig.Do.Subplot
                        figure;
                    end
                    
                    if isfield(FigConfig,'PlotFunction') && strcmpi(FigConfig.PlotFunction,'pcolor')
                        pcolor(f,1:narrays+1,derror1.'); shading flat;
                    else
                        surf(f,1:narrays+1,derror1.');  view(2); shading flat;
                    end
                    set(gca, 'XScale','log');
                    if FigConfig.Do.Title
                        title(['d_{error} in ' plane]);
                    end
                    ax = gca;
                    ax.XLabel.Interpreter = 'latex';
                    ax.YLabel.Interpreter = 'latex';
                    ax.TickLabelInterpreter = 'latex';
                    if FigConfig.Do.XLabel; xlabel('Frequency (Hz)'); end
                    ax.XLim = [500 min(f(end),fbandcutoff(end))];
                    ax.YLim =[1 narrays+1];
                    if FigConfig.Do.YTicks
                        ax.YTick = (1:narrays) + 1/2;
                        ax.YTickLabel = array_geo_disp;
                    else
                        ax.YTick = '';
                    end
                    ax.YDir = 'reverse';
                    if FigConfig.Do.XTicks
                        ax.XTick = foctaves;
                        ax.XTickLabel = foctavesTicks;
                    else
                        ax.XTick = '';
                    end
                    ax.CLim=[-30 0]; 
                    if FigConfig.Do.Cbar
                        cb=colorbar; 
                        set(get(cb,'YLabel'),'Interpreter','latex','string','$20\log_{10}d_{error}$ (dB)','fontsize',FigConfig.FontSize);
                        cb.TickLabelInterpreter = 'latex';
                    end
                    set(gca,'FontSize',FigConfig.FontSize);

        end
    end
end

   


end



