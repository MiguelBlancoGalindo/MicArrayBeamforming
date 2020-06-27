function [gamma,beta,r] = car2lat(x,y,z)
% Converts from cartesian to lateral angle coordinates

% Cartesian are defined as x to the right, y facing front and z up. 
% gamma, lateral angle in degrees from y clockwise. Takes values from -180 to 180
% beta, polar angle in degrees from intersection of the xy plane with 
% a parallel yz plane with the equivalent azimuth angle. Beta is the angle
% between the horizontal intersection and the point of interest at the
% parallel yz plane. Takes values from -90 to 90, positive facing upwards.
% r, radius of sphere from origin of coordinates

r = sqrt(x^2+y^2+z^2);
gamma = asin(x/r)*180/pi;
beta = atan2(z,y)*180/pi;

return;
