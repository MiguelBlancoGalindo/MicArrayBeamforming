function Config=checkConfigArrayMan(Config)
% function Config=checkConfigArrayMan(Config)
% function that checks that the configuration parameters of the array
% manifold are not missing or incorrect. 
%
% input and output arguments: 
%   Config: configuration structure containing all the settings.
%
% Note: The array manifold also depends on Config.Filter

if strcmp(Config.Array.Geo,'circular_rigid_sphere') || strcmp(Config.Array.Geo,'spherical_rigid_sphere')
    Config.ArrayMan.ManifoldType='spherical_scatterer'; %
elseif strcmp(Config.Array.Geo,'circular_rigid_cylinder') || strcmp(Config.Array.Geo,'stacked_circular_rigid_cylinder') ...
    || strcmp(Config.Array.Geo,'spiral')
    Config.ArrayMan.ManifoldType='cylindrical_scatterer'; %
else
    Config.ArrayMan.ManifoldType='freefield';
end


if ~isfield(Config.ArrayMan,'LookAng') || ~isfield(Config.ArrayMan.LookAng,'Az') || ~isfield(Config.ArrayMan.LookAng,'El')
if strcmp(Config.ArrayMan.SteerSpace,'2Daz')    
    % performer angle in degrees relative to array centre (2D only)
    SourceDirection = Config.ArrayMan.PerformerAngle;
    Config.ArrayMan.LookAng.Az = [-180+360/Config.ArrayMan.AzResolution:360/Config.ArrayMan.AzResolution:180].';
    Config.ArrayMan.LookAng.El = repmat(SourceDirection(1,1),Config.ArrayMan.AzResolution,1);
    naz = length(Config.ArrayMan.LookAng.Az);
	Config.ArrayMan.SolidAngle = 2*pi;
    Config.ArrayMan.Qweights = ones(naz,1)./naz;
elseif strcmp(Config.ArrayMan.SteerSpace,'2Del')
    % performer angle in degrees relative to array centre (2D only)
    SourceDirection = Config.ArrayMan.PerformerAngle;
    Config.ArrayMan.LookAng.El = [-180+360/Config.ArrayMan.ElResolution:360/Config.ArrayMan.ElResolution:180].';
    %Config.ArrayMan.LookAngEl = [-179:-91, -89:89, 91:180];
    Config.ArrayMan.LookAng.Az = repmat(SourceDirection(1,2),Config.ArrayMan.ElResolution,1);
    nel = length(Config.ArrayMan.LookAng.El);
    Config.ArrayMan.SolidAngle = 2*pi;
    Config.ArrayMan.Qweights = ones(nel,1)./nel;
elseif strcmp(Config.ArrayMan.SteerSpace,'3D')
    SA.SpkGrid.SamplingScheme = Config.ArrayMan.SteerGrid.Type;
    SA.SpkGrid.R = 1;
    if isfield(Config.ArrayMan.SteerGrid,'Points') && ~isempty(Config.ArrayMan.SteerGrid.Points)
        SA.SpkGrid.L = Config.ArrayMan.SteerGrid.Points;
    elseif isfield(Config.ArrayMan.SteerGrid,'Res') && ~isempty(Config.ArrayMan.SteerGrid.Res)
        SA.SpkGrid.L = 360*180/Config.ArrayMan.SteerGrid.Res.^2;
    else
        error('undefined number of point or resolution for sampling of the sphere');
    end
    SA = sphereSampling(SA,'Spk');
    
    Config.ArrayMan.LookAng.Az = 90-SA.SpkGrid.SteerDir(:,2);
    Config.ArrayMan.LookAng.Az(Config.ArrayMan.LookAng.Az>180) = Config.ArrayMan.LookAng.Az(Config.ArrayMan.LookAng.Az>180)-360;
    Config.ArrayMan.LookAng.Az(Config.ArrayMan.LookAng.Az<-179) = Config.ArrayMan.LookAng.Az(Config.ArrayMan.LookAng.Az<-179)+360;
    Config.ArrayMan.LookAng.El = 90-SA.SpkGrid.SteerDir(:,1);
    Config.ArrayMan.SolidAngle = 4*pi;
    Config.ArrayMan.Qweights = SA.SpkGrid.QWeights/sum(SA.SpkGrid.QWeights);
    
end
end

Config.ArrayMan.c=343; % assumed speed of sound
Config.Array.falias = Config.ArrayMan.c./(2*Config.Array.Spacing);

SNRIsodB = inf;     %Signal to isotropic noise ratio. Choose inf for no noise
DRRdB = inf;     %Direct-to-reverberant ratio. Choose inf for no diffuse noise
Config.ArrayMan.SNRIsodB = repmat(SNRIsodB,size(Config.Array.Offset));   %#ok<REPMAT>
Config.ArrayMan.DRRdB = repmat(DRRdB,size(Config.Array.Offset)); %#ok<REPMAT>


if ~isfield(Config.Filter,'FreqVector') || isempty(Config.Filter.FreqVector)
    error('Missing frequency vector');
end
if ~isfield(Config.Filter,'FreqMode') || isempty(Config.Filter.FreqMode)
    error('Missing frequency mode');
else
    if ~strcmpi(Config.Filter.FreqMode,'filter') && ~strcmp(Config.Filter.FreqMode,'discrete') 
        error(['Unknown FreqMode. Choose from ''', 'filter', ''' or ''' ,'discrete', '''']);
    end
end
if Config.Filter.FreqVector(end)*2>Config.Filter.Fs
    error('Frequency vector above Nyquist frequency')
end

end

