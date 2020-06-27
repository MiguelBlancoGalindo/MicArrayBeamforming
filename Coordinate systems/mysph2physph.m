function [phi_P,theta_P,r] = mysph2physph(phi_M,theta_M,r)
% Converts from my coordinate system to the standard physics coordinate
% system

% Cartesian are defined as x to the right, y facing front and z up. 
% r, radius of sphere from origin of coordinates
% My coordinate system
    % phi_M, azimuth in degrees from y clockwise. Takes values from -180 to 180
    % theta_M, elevation in degrees from xy plane upwards. Takes values from -90 to 90
    
% Standard physics coordinate system
    % phi_P, azimuth in degrees from x counterclockwise. 
    % theta_P, inclination angle in degrees from zenith downwards. 
    
% Test if wrapping beyond 360 is needed / recalculate values between -180
% to 180
phi_P = 90 - phi_M;
theta_P = 90 - theta_M;

end
