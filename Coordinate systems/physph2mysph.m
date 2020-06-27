function [phi_M,theta_M,r] = physph2mysph(phi_P,theta_P,r)
% Converts from the standard physics coordinate system to my coordinate
% system

% Cartesian are defined as x to the right, y facing front and z up. 
% r, radius of sphere from origin of coordinates

% Standard physics coordinate system
    % phi_P, azimuth in degrees from x counterclockwise. 
    % theta_P, inclination angle in degrees from zenith downwards. 
% My coordinate system
    % phi_M, azimuth in degrees from y clockwise. Takes values from -180 to 180
    % theta_M, elevation in degrees from xy plane upwards. Takes values from -90 to 90
    
% Test if wrapping beyond 360 is needed / recalculate values between -180
% to 180
phi_M = 90 - phi_P;
theta_M = 90 - theta_P;

end
