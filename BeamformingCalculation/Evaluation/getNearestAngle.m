function [el,az,ind,deltaOmega,deltaEl,deltaAz] = getNearestAngle(Config,theta0,phi0,varargin)
% function [el,az,ind,deltaOmega,deltaEl,deltaAz] = getNearestAngle(Config,theta0,phi0,varargin)
% function that calculates the nearest angle either in 2D or 3D from a
% specified angle with respect to the steering angle vector (default) or
% the microphone array positions, both included in Config. 
%
% input arguments:
%   Config: configuration struct containing all the settings.
%   theta0: elevation of the desired angle to find closest steering angle.
%   azimuth0: azimuth of the desired angle to find closest steering angle.
%
% optional input arguments:
%   'SPK': string to indicate that the closest angle to be found is with
%       respect to the steering angles (default is omitted). 
%   'MIC': string to indicate that the closest angle to be found is with
%       respect to the microphone array positions. 
%
% output arguments:
%   el: elevation angle of the closest angle.
%   az: azimuth angle of the closest angle.
%   ind: index of the (steering angle or microphone position) vector with 
%       closest angle.
%   deltaOmega: solid angle difference between target and closest angle.
%   deltaEl: difference in elevation between target and closest angle.
%   deltaAz: difference in azimuth between target and closest angle.


if nargin<4 || strcmpi(varargin{1},'SPK')
    phi = Config.ArrayMan.LookAng.Az;
    theta = Config.ArrayMan.LookAng.El;
    if length(phi)==1 
        phi = repmat(phi,length(theta),1);
    elseif length(theta)==1 
        theta = repmat(phi,length(theta),1);
    end
elseif strcmpi(varargin{1},'MIC')
    [phi,theta,~] = car2sph(Config.Array.MicPos(:,1),Config.Array.MicPos(:,2),Config.Array.MicPos(:,3));
end

L = length(theta);
dist = nan(L,1);
distaz = nan(L,1);
distel = nan(L,1);

for l=1:L
    [dist(l),distaz(l),distel(l)] = spherror(phi0,theta0,phi(l),theta(l));
end

[deltaOmega,ind] = min(dist);
deltaAz = distaz(ind);
deltaEl = distel(ind);

az = phi(ind);
el = theta(ind);