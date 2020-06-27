function [gamma,beta,r] = sph2lat(phi,theta,r)
% Converts from spherical to lateral angle coordinates

% phi, azimuth in degrees from -180 to 180
% theta, elevation in degrees from -90 to 90
% gamma, lateral angle in degrees from from -180 to 180
% beta, polar angle in degrees from -90 to 90

[x,y,z] = sph2car(phi,theta,r);
[gamma,beta] = car2lat(x,y,z);

return;
