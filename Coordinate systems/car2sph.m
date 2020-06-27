function [phi,theta,r] = car2sph(x,y,z)
% Converts from cartesian to spherical coordinates

% Cartesian are defined as x to the right, y facing front and z up. 
% r, radius of sphere from origin of coordinates
% phi, azimuth in degrees from y clockwise. Takes values from -180 to 180
% theta, elevation in degrees from xy plane upwards. Takes values from -90 to 90

r = sqrt(x.^2 + y.^2 + z.^2);
phi = atan2(x,y).*180/pi;
theta = asin(z./r).*180/pi;

return;% theta, phi;
