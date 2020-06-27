function Config = configArray(ArrayGeo,ArrayGeoDisp,ArrayCriterion,M,r,SensorOffset,SensorOffsetStr,SensorSNR,SensorNoiseStr,CalErrorMagdB,CalErrorFreqdB, varargin)
% function Config = configArray(ArrayGeo,ArrayGeoDisp,ArrayCriteria,M,r,SensorOffset,SensorOffsetStr,SensorSNR,SensorNoiseStr,CalErrorMagdB,CalErrorFreqdB, varargin)
% function that configures the microphone array based on the properties
% passed as input parameters.
% input parameters:
%   ArrayGeo: geometry of microphone array. Choose from (see "array_setup.m"): 
%       'linear_bs' (broadside), 'linear_ef' (endfire), 'circular', 
%       'dualcircular', 'multicircular', 'square', 'rectangular', 
%       'circular_grid', 'circular_rigid_sphere', 'circular_rigid_cylinder',
%       'stacked_circular' (vertically stacked circular), 
%       'stacked_circular_rigid_cylinder', 'spherical' or 
%       'spherical_rigid_sphere'.
%   ArrayGeoDisp: sortened version of ArrayGeo for display purposes only.
%   ArrayCriterion: criterion to compare different geometries. Choose from
%       'aperture' to fix the aperture and number of microphones, 'spacing'
%       to fix the spacing and microphones or 'aperture_and_spacing' to fix
%       the aperture and spacing but not the number of microphones.
%   M: number of microphones.
%   r: radius of equivalent maximum aperture of 2r.
%   SensorOffset: magnitude of microphone positioning offset.
%   SensorOffsetStr: string describing microphone positoning offset. 
%   SensorSNR: magnitude of SNR (in dB) from which sensor noise is derived.
%   SensorNoiseStr: string describing microphone self noise. 
%   CalErrorMagdB: magnitude of calibration error in dB.
%   CalErrorFreqdB: magnitude calibration deviation as a function of 
%        frequency (in dB).
%
% output arguments: 
%   Config: configuration structure containing all the settings.

Config.Array.Geo = ArrayGeo;
Config.Array.GeoDisp = ArrayGeoDisp;
Config.Array.Criteria = ArrayCriterion;
Config.Array.M = M;
Config.Array.r = r;

%Errors in Array Manifold. If more than one type of error is to be
%evaluated they need to be included with all combinatons to be assessed in
%a vector of the same length for each of the error types
Config.Array.Offset = SensorOffset;
Config.Array.OffsetStr = SensorOffsetStr;
Config.Array.SNR = SensorSNR;
Config.Array.NoiseStr = SensorNoiseStr;
Config.Array.CalErrorMagdB = CalErrorMagdB;
Config.Array.CalErrorFreqdB = CalErrorFreqdB;

iinputMicPos = 0;
doArrayMan = false;
nVarargs = length(varargin);
for n=1:2:nVarargs
    if strcmpi(varargin{n},'MicPos')
        MicPos = varargin{n+1};
        iinputMicPos = 1;
    elseif strcmpi(varargin{n},'NCircles')
        Config.Array.NCircles = varargin{n+1};
    elseif strcmpi(varargin{n},'StackedCircleSpacing')
        Config.Array.StackedCircleSpacing = varargin{n+1};
    elseif strcmpi(varargin{n},'MultiCircleAlignment')
        Config.Array.MultiCircleAlignment = varargin{n+1};
    elseif strcmpi(varargin{n},'SourceDirection')
        doArrayMan = true;
        SourceDirection = varargin{n+1};
    end
end
if iinputMicPos == 0
    Config = initialise_array_setup(Config);
    Config = array_setup(Config);
else
    Config.Array.MicPos = MicPos;
end

Config.Array.Spacing = calculate_array_spacing('ArrayGeo',Config.Array.Geo, 'ArrayPos', Config.Array.MicPos);
Config.Array.Aperture = array_aperture(Config.Array.MicPos, Config.Array.Geo, 'far_field');

if doArrayMan
Config=configArrayMan(Config,SourceDirection);
end

end
