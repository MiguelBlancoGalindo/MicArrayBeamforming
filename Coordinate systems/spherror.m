function [errarc,az_err,el_err] = spherror(u_phi,u_theta,v_phi,v_theta)
% Calculates, in degrees, the arc length (errarc) between two 2D directions
% and resolves it into horizontal (err1) and vertical (err2) components
% according to either spherical or lateral angle coordinates
% u_phi, azimuth angle of reference point
% u_theta, elevation angle of reference point
% u_phi, azimuth angle of estimated point
% u_theta, elevation angle of estimated point


% Conversion of u and v to cartesian coordinates
u = zeros(size(u_phi,1),size(u_phi,2));
v = zeros(size(v_phi,1),size(v_phi,2));
[u(:,1), u(:,2), u(:,3)] = sph2car(u_phi,u_theta,1);
[v(:,1), v(:,2), v(:,3)] = sph2car(v_phi,v_theta,1);

%projection of v over u in azimuth
w = zeros(size(u_phi,1),size(u_phi,2));
[w(:,1), w(:,2), w(:,3)] = sph2car(v_phi,u_theta,1);
%vectors from centre of circumference at z equal to that of vectors u and
%w
u1 = u - [zeros(size(u,1),1), zeros(size(u,1),1), u(:,3)];
w1 = w - [zeros(size(u,1),1), zeros(size(u,1),1), w(:,3)];
%projection of v over u in elevation
x = zeros(size(u,1),size(u,2));
[x(:,1), x(:,2), x(:,3)] = sph2car(u_phi,v_theta,1);
  
% Calculate the greater circle arc length
errarc = acos(dot(u,v,2));

% Find the error components
if  errarc > 0
    u1norm = sqrt(u1(:,1).^2 + u1(:,2).^2 + u1(:,3).^2);
    w1norm = sqrt(w1(:,1).^2 + w1(:,2).^2 + w1(:,3).^2);
    az_err = abs(u1norm.*acos(dot(u1,w1,2)./(u1norm.*w1norm)));
    el_err = abs(acos(dot(u,x,2))); 

else
    az_err = 0;
    el_err = 0;
end

% Express angular errors in degrees
errarc = errarc*180/pi;
az_err = az_err*180/pi;
el_err = el_err*180/pi;

end
