% script to evaluate the performance of different arrays and beamformers in
% 3D.

addpath(genpath('../'));

ArrayGeo={'circular_rigid_cylinder','spherical_rigid_sphere'};
ArrayGeoDisp={'C-RC','S-RS'};

%number of microphones
M=32;
%radius of the equivalent circular/spherical array with the same aperture/spacing depending on array_geo_comparsison
r=0.1;
rmm = r*1000;   %radius in mm
% ArrayCriteria choose from 'aperture', 'spacing' or 'aperture_and_spacing'
ArrayCriteria = 'aperture';

%SensorOffset choose from 0, 0.001 or 0.002 or [] for ideal case regardless
%of other parameter values
SensorOffset = [];  
%SensorOffsetStr choose from 'no', '1mm' or '2mm'
if isempty(SensorOffset) || SensorOffset==0
    SensorOffsetStr = 'no';
else
    SensorOffsetStr = [num2str(SensorOffset) 'mm'];
end
%SensorSNR. Choose inf for no noise
SensorSNR = inf;

if SensorSNR==inf
    SensorNoiseStr = 'no';
else
    SensorNoiseStr = ['-' num2str(SensorSNR) 'dB'];
end

%frequency independent calibration error magnitude in dB (0 for no error)
CalErrorMagdB = 0;
%frequency independent calibration error magnitude in dB (0 for no error)
CalErrorFreqdB = 0;

% Source and target directions and amplitudes
SourceDirection=[0 0];    % nSourceDir x 2 where each row is el and az for each target
Interferers=[0 60];   % nSourceDir x (2*nInt) where each row are pairs of el and az interferers next to each other

narrays = length(ArrayGeo);
%% configuring microphone array
ConfigAllArrays = cell(narrays,1);

for iarray=1:narrays
    %use this to create MicPos
    ConfigAllArrays{iarray} = configArray(ArrayGeo{iarray},ArrayGeoDisp{iarray},ArrayCriteria,M,r,SensorOffset,SensorOffsetStr,SensorSNR,SensorNoiseStr,CalErrorMagdB,CalErrorFreqdB,SourceDirection);
    %use this if already have MicPos
    %ConfigAllArrays{iarray,iM,ir} = configArray(ArrayGeo{iarray},ArrayGeoDisp{iarray},ArrayCriteria,M(iM),r(ir),SensorOffset,SensorOffsetStr,SensorSNR,SensorNoiseStr,CalErrorMagdB,CalErrorFreqdB,SourceDirection,'MicPos',MicPos{iarray,iM,ir});
end
      
%% configuring array manifold

Config = ConfigAllArrays{1};
TargetDirection = [0 0];
Config.ArrayMan.LookDistance='farfield'; % 'farfield', 'focalpoint'
Config.ArrayMan.SteerSpace='3D'; % '2Daz', '2Del', '3D';
Config.ArrayMan.SteerGrid.Type = 'Lebedev'; %'Lebedev', 'Equiangular'
Config.ArrayMan.SteerGrid.Points = 590; %590, 360/5*180/5
%Config.ArrayMan.AzResolution=360; % Resolution in azimuth (sample points around the circle)
%Config.ArrayMan.ElResolution=1;   % Resolution in elevation 
Config.ArrayMan.PerformerDistance=1; % distance for 'focalpoint' array manifold
Config.ArrayMan.PerformerAngle=TargetDirection;

Config.Filter.Fs=48000; % filter fs
Config.Filter.FreqMode='discrete'; % 'filter' (narrowband processing for each FFT bin; 'discrete' (use specified vector)

if strcmp(Config.Filter.FreqMode,'filter')
    Config.Filter.Nfft = 1024;
    NfftNyquist = floor(Config.Filter.Nfft/2)+1;
    Config.Filter.FreqVector = 0:Config.Filter.Fs/Config.Filter.Nfft:Config.Filter.Fs/Config.Filter.Nfft*(NfftNyquist-1);
elseif strcmp(Config.Filter.FreqMode,'discrete')
    %Config.Filter.FreqVector=20:20:20000;
    nFreqsPerOctave = 20;
    Config.Filter.FreqVector = 1000*2.^(-5-1/2:1/nFreqsPerOctave:4+1/2);
end

for iarray=1:narrays
    ConfigAllArrays{iarray}.ArrayMan = Config.ArrayMan;
    ConfigAllArrays{iarray}.Filter = Config.Filter;
end

%% calculating array manifold
ArrayManAllArrays = cell(narrays,1);
for iarray=1:narrays
    disp(['Microphone array: ' ConfigAllArrays{iarray}.Array.GeoDisp]);
    [ArrayMan,Config]=GetArrayManifold(ConfigAllArrays{iarray}); 
    ConfigAllArrays{iarray} = Config;
    ArrayManAllArrays{iarray} = ArrayMan;
end

%% configuring beamforming calculation and evaluation
TargetDirection = [0 0];    %elevation and azimuth
InterfererDirection = [];   %elevation and azimuth
Config = ConfigAllArrays{1};
% choose beamformer methods from 'ds','sda','mvdr','lcmv','acc','ls'
Config.Filter.Method={'ls'};
N = [3 5];
Config.Filter.ls.mode='hypercardioid';
Config.Filter.ls.order=N;
Config.Filter.ls.EQ = false;

Config.Filter.RegMethod='wnglimit';
Config.Filter.RegParam=-10;% -10 dB WNGmin;

Config.Scene.Targets.Gain = 1;  
Config.Scene.Interferers.Gain = 1; 
Config.Scene.Diffuse.Coeff = 1;
Config.Scene.Targets.Angle=TargetDirection;
Config.Scene.Interferers.Angle=InterfererDirection; % angles of interfering sources, column for multiple interferers
Config.Scene.AngRange = [0 60];
SetupAudio = [];
Rxx = []; 
Config=checkConfigFilters(Config);
narrays = size(ConfigAllArrays,1);
for iarray = 1:narrays
    ConfigAllArrays{iarray}.Filter = Config.Filter;
    ConfigAllArrays{iarray}.Scene = Config.Scene;
end

%% Beamforming calculation and evaluation
[BFEvalAllArrays, BFWeightsAllArrays, ConfigAllArrays] = calculate_beamforming_weights_and_evaluation_different_arrays(ConfigAllArrays,ArrayManAllArrays);

%% rearrange struct with new format
BFEvalAllArrays2 = rearrange_directional_response_measures(ConfigAllArrays,BFEvalAllArrays);

%% ploting directional responses in 3D
FigConfig.Do.Save=0;
FigConfig.Arrays = {'all'};
FigConfig.Methods = {'all'};
FigConfig.Freq = [125,250,500,1000,2000,4000,8000,16000];
FigConfig.NoiseType = {''};
FigConfig.Noises = 1;
FigConfig.Do.OnAxis.Phase = true;
FigConfig.Do.OnAxis.Mag = true;
FigConfig.Do.Az = true;
FigConfig.Do.El = false;
FigConfig.Do.Title = false;
FigConfig.Do.Legend = false;
FigConfig.Path = [];
PlotDirectionalResponse3DAllArrays(ConfigAllArrays,BFEvalAllArrays2, FigConfig);

%% plotting measures of directional responses for all arrays
FigConfig.Do.Save=0;
FigConfig.Do.NewFig = true;
FigConfig.Do.Subplot = false;
FigConfig.Do.XLabel = true; 
FigConfig.Do.YLabel = true;
FigConfig.Do.XTicks = true; 
FigConfig.Do.YTicks = true;
FigConfig.Compare = 'methods';   %plot 'arrays' or 'methods' in same plot
%FigConfig.Arrays = {'all'};
FigConfig.Arrays = {'circular'};
%FigConfig.Arrays = {'circular_rigid_cylinder'};
%FigConfig.Arrays = {'circular','dualcircular','stacked_circular','spherical','circular_rigid_cylinder','stacked_circular_rigid_cylinder','circular_rigid_sphere','spherical_rigid_sphere'};
FigConfig.Methods = {'all'};
%FigConfig.Methods = {'sda'};
FigConfig.Metrics = {'DI'};    %'BW','SSL','DI','WNG'
FigConfig.FRange = 'fullband';     %'fullband','operative','nonaliased' 
FigConfig.NoiseType = {'-10dBWNG_360el'};
FigConfig.Type = 'absolute';  %difference or absolute
FigConfig.Do.Title = false;
FigConfig.Do.Legend = false;
FigConfig.Do.Zoom = false;
%plot_directional_response_measures_all_arrays(ConfigAllArrays,BFEvalABFEllArrays, FigConfig);
PlotDirectionalResponseMeasuresAllArrays(ConfigAllArrays,BFEvalAllArrays2, FigConfig);
%% plotting frequency range for all arrays and beamformers to choose
FigConfig.Do.Save=0;
FigConfig.Arrays = {'all'};
FigConfig.Methods = {'ds','sda','ls'};
PlotFrequencyRange(BFEvalAllArrays2,ConfigAllArrays,fRangeLSBAz, FigConfig)
