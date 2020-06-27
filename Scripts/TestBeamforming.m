%%%%%%%%%%%%%%%%%%% TestBeamforming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a script that serves as an example of the steps for the
% beamforming toolbox, from the microphone array design, array manifold
% calculation, beamforming computation and analysis of the performance
% metrics
addpath(genpath('../'));

%%% Defining the array geometries to be evaluated
% Choose from 'linear_ef', 'linear_bs', 'rectangular', 'circular','dualcircular','spherical',...
%    'circular_rigid_cylinder','circular_rigid_sphere','spherical_rigid_sphere'
ArrayGeo={'linear_ef', 'circular'};
% equivalent abbreviations (can have any name, for display purposes only)
ArrayGeoDisp={'LE','C'};

% Settings for all arrays
M=32;   % number of microphones
r=0.1;  % radius of equivalent maximum aperture

% ArrayCriteria choose from 'aperture', 'spacing' or 'aperture_and_spacing'
% In this case maximum aperture of all arrays is fixed. 
ArrayCriteria = 'aperture';

% Sensor positioning offset. Choose a value (in m). If [] is used, ideal
% array manifold is used regardless of other parameter values in SensorSNR,
%  CalErrorMagdB or CalErrorFreqdB (overridden).
SensorOffset = [];  
%String that identifies the value included in SensorOffset.
if isempty(SensorOffset) || SensorOffset==0
    SensorOffsetStr = 'no';
else
    SensorOffsetStr = [num2str(SensorOffset) 'mm'];
end
% SensorSNR. Choose inf for no noise
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

narrays = length(ArrayGeo);

%% configuring microphone array
ConfigAllArrays = cell(narrays,1);

for iarray=1:narrays
    %use this to create MicPos from ArrayGeo
    ConfigAllArrays{iarray} = configArray(ArrayGeo{iarray},ArrayGeoDisp{iarray},ArrayCriteria,M,r,SensorOffset,SensorOffsetStr,SensorSNR,SensorNoiseStr,CalErrorMagdB,CalErrorFreqdB);
    %use this if already have MicPos
    %ConfigAllArrays{iarray} = configArray(ArrayGeo{iarray},ArrayGeoDisp{iarray},ArrayCriteria,M,r,SensorOffset,SensorOffsetStr,SensorSNR,SensorNoiseStr,CalErrorMagdB,CalErrorFreqdB,'MicPos',MicPos{iarray,iM,ir});
end
ConfigArray = ConfigAllArrays{1};

%% configuring array manifold

Config = ConfigArray;
TargetDirection = [0 0];    %[Elevation azimuth] 
Config.ArrayMan.LookDistance='farfield'; % 'farfield', 'focalpoint'
Config.ArrayMan.SteerSpace='2Daz'; % '2Daz', '2Del', '3D';
% Config.ArrayMan.SteerGrid.Type = 'Lebedev'; %'Lebedev', 'Equiangular', 'Gaussian', 'Uniform'
% Config.ArrayMan.SteerGrid.Points = 590; %590, 360/5*180/5
Config.ArrayMan.AzResolution=360; % Resolution in azimuth (sample points around the circle)
Config.ArrayMan.ElResolution=1;   % Resolution in elevation 
Config.ArrayMan.PerformerDistance=1; % distance for 'focalpoint' array manifold
Config.ArrayMan.PerformerAngle=TargetDirection;

Config.Filter.Fs=48000; % filter fs
Config.Filter.FreqMode='discrete'; % 'filter' (narrowband processing for each FFT bin; 'discrete' (use specified vector)

if strcmp(Config.Filter.FreqMode,'filter')
    Config.Filter.Nfft = 1024;
    NfftNyquist = floor(Config.Filter.Nfft/2)+1;
    Config.Filter.FreqVector = 0:Config.Filter.Fs/Config.Filter.Nfft:Config.Filter.Fs/Config.Filter.Nfft*(NfftNyquist-1);
elseif strcmp(Config.Filter.FreqMode,'discrete')
    Config.Filter.FreqVector=20:20:20000;
    %nFreqsPerOctave = 20;
    %Config.Filter.FreqVector = 1000*2.^(-5-1/2:1/nFreqsPerOctave:4+1/2);
end

for iarray=1:narrays
    ConfigAllArrays{iarray}.ArrayMan = Config.ArrayMan;
    ConfigAllArrays{iarray}.Filter = Config.Filter;
end

%% calculating array manifold steering matrix

ArrayManAllArrays = cell(narrays,1);
for iarray=1:narrays
    disp(['Microphone array: ' ConfigAllArrays{iarray}.Array.GeoDisp]);
    [ArrayMan,Config]=GetArrayManifold(ConfigAllArrays{iarray}); 
    ConfigAllArrays{iarray} = Config;
    ArrayManAllArrays{iarray} = ArrayMan;
end

%% configuring beamforming calculation and evaluation
TargetDirection = [0 0];    %elevation and azimuth
InterfererDirection = [0 60];   %elevation and azimuth
Config = ConfigAllArrays{1};
% choose beamformer methods from 'ds','sda','mvdr','lcmv','acc','ls'
Config.Filter.Method={'ds','sda','mvdr','ls'};
N = 4;
Config.Filter.ls.mode='hypercardioid';  %'hypercardioid', 'cardioid', 'chebyshev'
Config.Filter.ls.order=N;   
Config.Filter.ls.EQ = false;

Config.Filter.RegMethod='wnglimit';     %'external' or 'wnglimit'
Config.Filter.RegParam=-10;             % -10 dB WNGmin;

Config.Scene.Targets.Gain = 1;  
Config.Scene.Interferers.Gain = 1; 
Config.Scene.Diffuse.Coeff = 1;
Config.Scene.Targets.Angle=TargetDirection; %elevation and azimuth
Config.Scene.Interferers.Angle=InterfererDirection; %elevation and azimuth
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

%% ploting directional responses in azimuth and elevation
FigConfig.Do.Save=0;
FigConfig.Arrays = {'all'};
FigConfig.Methods = {'all'};
FigConfig.NoiseType = {''};
FigConfig.Noises = 1;
FigConfig.Do.Title = true;
FigConfig.Do.Legend = false;
%FigConfig.Path = '/user/HS228/mb00859/Dropbox/Surrey Uni/PhD/Documents/Journal articles/Draft/';
FigConfig.Path = [];
PlotDirectionalResponse2DAllArrays(ConfigAllArrays,BFEvalAllArrays2, FigConfig);
%% plotting measures of directional responses for all arrays
FigConfig.Do.Save=0;
FigConfig.Do.NewFig = true;
FigConfig.Do.Subplot = false;
FigConfig.Do.XLabel = true; 
FigConfig.Do.YLabel = true;
FigConfig.Do.XTicks = true; 
FigConfig.Do.YTicks = true;
FigConfig.Compare = 'arrays';   %plot 'arrays' or 'methods' in same plot
FigConfig.Arrays = {'all'};
%FigConfig.Arrays = {'circular_rigid_cylinder'};
%FigConfig.Arrays = {'circular','dualcircular','stacked_circular','spherical','circular_rigid_cylinder','stacked_circular_rigid_cylinder','circular_rigid_sphere','spherical_rigid_sphere'};
%FigConfig.Methods = {'all'};
FigConfig.Methods = {'sda','mvdr'};
FigConfig.Metrics = {'BW'};    %'BW','SSL','DI','WNG'
FigConfig.FRange = 'fullband';     %'fullband','operative','nonaliased' 
FigConfig.NoiseType = {'-10dBWNG_360el'};
FigConfig.Type = 'absolute';  %difference or absolute
FigConfig.Do.Title = true;
FigConfig.Do.Legend = true;
FigConfig.Do.Zoom = false;
FigConfig.Do.NewFig = true;
%plot_directional_response_measures_all_arrays(ConfigAllArrays,BFEvalABFEllArrays, FigConfig);
PlotDirectionalResponseMeasuresAllArrays(ConfigAllArrays,BFEvalAllArrays2, FigConfig);
%% plotting frequency range for all arrays and beamformers to choose
% calculating the directional response error between synthesised and target 
% responses for least-squares beamformer (LSB).
[derror,~,fRangeLSB] = directionalResponseErrorLSB(BFEvalAllArrays2, ConfigAllArrays, 0);

clear FigConfig
FigConfig.Do.Save=0;
FigConfig.Do.Subplot = false;
FigConfig.Do.Title = false;
FigConfig.Do.XLabel = true; 
FigConfig.Do.YTicks = true;
FigConfig.Do.XTicks = true; 
FigConfig.Do.Legend = true;
FigConfig.Arrays = {'all'};
FigConfig.Methods = {'ds','sda','ls'};
figure; 
PlotFrequencyRange(BFEvalAllArrays2,ConfigAllArrays,fRangeLSB, FigConfig);

%% plotting directional response error between synthesised and target responses for LSB

FigConfig.Do.Cbar = true;
figure; PlotLSBDirRespError(BFEvalAllArrays2,ConfigAllArrays, derror, FigConfig);

