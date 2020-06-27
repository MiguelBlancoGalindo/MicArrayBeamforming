function [x,y,z] = cyl2car(phi,z,rho)
% Converts from cylindrical coordinates to cartesian

% Cartesian are defined as x to the right, y facing front and z up. 
% rho, radius of cylinder from origin of coordinates
% phi, azimuth in degrees from y clockwise. Takes values from -180 to 180
% z, height of cylinder

x = rho .* sin(phi*pi/180);
y = rho .* cos(phi*pi/180);
z = z;

return;% x, y, z;
