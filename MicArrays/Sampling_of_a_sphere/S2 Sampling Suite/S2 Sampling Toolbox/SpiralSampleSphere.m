function [V,Tri]=SpiralSampleSphere(N,vis)
% Produce an approximately uniform sampling of a unit sphere using a spiral
% method [1]. According to this approach particle (i.e., sample) longitudes 
% are proportional to particle rank (1 to N) and latitudes are assigned to
% ensure uniform sampling density.  
%
% INPUT:
%   - N     : desired number of particles. N=200 is the default setting.
%             Note that sampling becomes more uniform with increasing N. 
%   - vis   : optional logical input argument specifying if you the 
%             spiral-based sampling should be visualized. vis=false is the 
%             default setting.
%
% OUTPUT:  
%   - V     : N-by-3 array of vertex (i.e., sample) co-ordinates.
%   - Tri   : M-by-3 list of face-vertex connectivities. 
%
%
% REFERENCES:
% [1] Christopher Carlson, 'How I Made Wine Glasses from Sunflowers', 
%     July 8, 2011. url: http://blog.wolfram.com/2011/07/28/how-i-made-wine-glasses-from-sunflowers/
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: March, 2015  
%

if nargin<1 || isempty(N), N=200; end
if nargin<2 || isempty(vis), vis=false; end

N=round(N(1));

gr=(1+sqrt(5))/2;       % golden ratio
ga=2*pi*(1-1/gr);       % golden angle

i=0:(N-1);              % particle (i.e., point sample) index
lat=acos(1-2*i/(N-1));  % latitude is defined so that particle index is proportional to surface area between 0 and lat
lon=i*ga;               % position particles at even intervals along longitude

% Convert from spherical to Cartesian co-ordinates
x=sin(lat).*cos(lon);
y=sin(lat).*sin(lon);
z=cos(lat);
V=[x(:) y(:) z(:)];

% Is triangulation required?
Tri=[];
if nargout>1, Tri=fliplr(convhulln(V)); end

% Visualize the result
if ~vis, return; end
tr=SubdivideSphericalMesh(IcosahedronMesh,5);
figure('color','w')

if ~isempty(Tri), ha1=subplot(1,2,1); end
h=patch('faces',tr.Triangulation,'vertices',tr.X);
set(h,'EdgeColor','none','FaceColor',[0 0.8 0],'SpecularStrength',0.5)
axis equal off vis3d
hold on

n=2;
col=[1 0 0; 0 0 0];
for i=1:n
    plot3(x(i:n:N),y(i:n:N),z(i:n:N),'.r','MarkerSize',max(min(30*sqrt(1E3/N),30),5),'MarkerEdgeColor',col(i,:))
end

if N<300
    r=1.001;
    i=linspace(0,N-1,1E5);
    lat=acos(1-2*i/(N-1));
    lon=i*ga;
    x=sin(lat).*cos(lon);
    y=sin(lat).*sin(lon);
    z=cos(lat);
    plot3(r*x,r*y,r*z,'-k','Color',[0 0 0.4])
end

light 
lighting phong
view([20 30])

if isempty(Tri), return; end
ha2=subplot(1,2,2);

h=patch('faces',Tri,'vertices',V);
set(h,'FaceColor','w','EdgeColor','k')
axis equal off vis3d
hold on
for i=1:n
    plot3(V(i:n:N,1),V(i:n:N,2),V(i:n:N,3),'.r','MarkerSize',max(min(30*sqrt(1E3/N),30),5),'MarkerEdgeColor',col(i,:))
end
set(ha2,'CameraViewAngle',get(ha1,'CameraViewAngle'))
view([20 30])

