function [x,y,z] = sph2car(phi,theta,r)
% Converts from spherical coordinates to cartesian

% Cartesian are defined as x to the right, y facing front and z up. 
% r, radius of sphere from origin of coordinates
% phi, azimuth in degrees from y clockwise. Takes values from -180 to 180
% theta, elevation in degrees from xy plane upwards. Takes values from -90 to 90

x = r .* cos(theta*pi/180) .* sin(phi*pi/180);
y = r .* cos(theta*pi/180) .* cos(phi*pi/180);
z = r .* sin(theta*pi/180);

return;% x, y, z;
