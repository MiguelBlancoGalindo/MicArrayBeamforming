function [phi,theta,r] = lat2sph(gamma,beta,r)
% Converts from lateral angle to spherical coordinates

% Cartesian are defined as x to the right, y facing front and z up. 
% gamma, lateral angle in degrees from y clockwise. Takes values from -180 to 180
% beta, polar angle in degrees from intersection of the xy plane with 
% a parallel yz plane with the equivalent azimuth angle. Beta is the angle
% between the horizontal intersection and the point of interest at the
% parallel yz plane. Takes values from -90 to 90, positive facing upwards.
% r, radius of sphere from origin of coordinates
% phi, azimuth in degrees from y clockwise. Takes values from -180 to 180
% theta, elevation in degrees from xy plane upwards. Takes values from -90 to 90

[x,y,z] = lat2car(gamma,beta,r);
[phi,theta,~] = car2sph(x,y,z);

return;
