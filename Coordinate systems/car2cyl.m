function [phi,z,rho] = car2cyl(x,y,z)
% Converts from cartesian to cylindrical coordinates

% Cartesian are defined as x to the right, y facing front and z up. 
% rho, radius of cylinder from origin of coordinates
% phi, azimuth in degrees from y clockwise. Takes values from -180 to 180
% z, height of cylinder

rho = sqrt(x.^2 + y.^2);
phi = atan2(x,y).*180/pi;
z = z;

return;% theta, phi;
