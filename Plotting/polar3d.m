function polar3d(d,Omega,FigConfig)
% function polar3d(d,Omega,FigConfig)
% function to create a 3D polar plot of a given directional response.
% Adapted from Lorenzo Chiesi.
%
% input arguments:
%   d: directional response in 3d stacked in column vector (for single
%       frequency).
%   Omega: steering angle matrix (SurfacePoint x [theta,phi]).
%   FigConfig: struct to pass plotting parameters.


if size(d,2)>size(d,1) && size(d,1)==1
    d = d.';
end

if (size(d,1) < 500)
    %Create Theta/Phi mesh with 5deg resolution
    [phi,theta] = meshgrid( -180:5:180 , -90:5:90);
else
    %Create Theta/Phi mesh with 1deg resolution
    [phi,theta] = meshgrid( -180:180 , -90:90);
end

%Complete matrix wrapping around on phi dimension
Omega2 = Omega;
Omega2(:,2) = Omega2(:,2)+360;
Omega = [Omega; Omega2];
d = [d; d];

%Interpolate sparse point finding magnitude and phase value on grid
if isfield(FigConfig.Do,'dB') && FigConfig.Do.dB==true 
    Magnitude = griddata(Omega(:,2),Omega(:,1),db(d),phi,theta,'cubic');
else
    Magnitude = griddata(Omega(:,2),Omega(:,1),abs(d),phi,theta,'cubic');
end
Phase = griddata(Omega(:,2),Omega(:,1),abs(180/pi*angle(d)),phi,theta);
% Magnitude = griddata(Omega(:,1),Omega(:,2),abs(d),theta,phi,'cubic');
% Phase = griddata(Omega(:,1),Omega(:,2),abs(angle(d)),theta,phi);

%Convert to cartesian
[X,Y,Z] = sph2car(phi,theta,Magnitude); 

surf(X,Y,Z,Phase); shading flat;
colormap(gca,flipdim(jet,1));
caxis([0,180]);

if isfield(FigConfig,'Lim') && ~isempty(FigConfig.Lim)
    dmin = FigConfig.Lim(1);
    dmax = FigConfig.Lim(2);
else
    dmin = -max(abs(d));
    dmax = max(abs(d));
end

xlim([dmin, dmax]);
ylim([dmin, dmax]);
zlim([dmin, dmax]);
ax = gca;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';

axis vis3d

end
