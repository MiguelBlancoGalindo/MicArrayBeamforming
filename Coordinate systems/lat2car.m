function [x,y,z] = lat2car(gamma,beta,r)
% Converts from lateral angle coordinates to cartesian

% Cartesian are defined as x to the right, y facing front and z up. 
% gamma, lateral angle in degrees from y clockwise. Takes values from -180 to 180
% beta, polar angle in degrees from intersection of the xy plane with 
% a parallel yz plane with the equivalent azimuth angle. Beta is the angle
% between the horizontal intersection and the point of interest at the
% parallel yz plane. Takes values from -90 to 90, positive facing upwards.
% r, radius of sphere from origin of coordinates

x = r*sin(gamma*pi/180); 
y = r*cos(gamma*pi/180)*cos(beta*pi/180);
z = r*cos(gamma*pi/180)*sin(beta*pi/180);

return;