% coversphere3.m
%
% evenly distributes points on a circle on the surface of a sphere
%
% 160819PJ: based on coversphere2, developed to plot detailed geometry


% Switches
isCoords = 1;

% Parameters
r1 = 0.083;% 17cm sphere internal diameter
r2 = 0.085;% 17cm sphere external diameter
mRad = 0.5 * 0.0046;% Countryman B3 radius
gRad = 0.5 * 0.0102;% RS 666-4596 grommet radius
cRad1 = 0.5 * 0.0208;% RS 666-4622 grommet inside radius
cRad2 = 0.5 * 0.0300;% RS 666-4622 grommet outside radius
fRad1 = 0.5 * 0.006;% fixing inside radius
fRad2 = 0.5 * 0.010;% fixing outside radius
iTh = [0:5:360];
aColour = [0.8 0.2 0.2];% hemisphere join colour
rColourx = [0.0 0.0 1.0 ];% microphone colour
rColourg = [0.5 0.5 0.95 ];% grommet colour
rColourf = [0.5 0.95 0.5 ];% fixing colour
hAz = -0.5 * 360 / 32;% -5.625 deg, half of 11.25 deg
hEl = 0.0 * 360 / 24;% 7.5 deg, half of 15 deg
cEl = -60;% cable duct at the rear
% cAz = 180;% cable duct at the rear, just beside the microphone stand
cAz = 180 - (360 * 0.026 / (2 * pi * r2 * cos(cEl*pi/180)));
fAz = 180 * [1; 1];
fEl = [(12.6 - 60.0); (-13.2 - 60.0)];
nF = 2;

a = load('m32circle.txt');
[nM,nC] = size(a);
mAz = a(:,2);
mEl = a(:,3);
x = r2 * cos(mEl*pi/180) .* cos(mAz*pi/180);
y = r2 * cos(mEl*pi/180) .* sin(mAz*pi/180);
z = r2 * sin(mEl*pi/180);


figure(1);
cla;
axis equal;

% Plot sphere
[xs,ys,zs] = sphere(72);
hs = surf(0.999*r2*xs,0.999*r2*ys,0.999*r2*zs);
set(hs,'Facecolor',[0.95 0.95 0.75]);
l = light;
lighting phong;
hold on;

% Plot microphones
for iM = 1:nM,
  plot3(x(iM),y(iM),z(iM),'r.');
  aUnit = [x(iM); y(iM); z(iM)] / r2;% radial unit vector
%  if (abs(mEl(iM)) == 90),
%    bUnit = [0; 1; 0];
%  else,
    bUnit = [-sin(mAz(iM)*pi/180); cos(mAz(iM)*pi/180); 0];
%  end;% if
  cUnit = cross(aUnit, bUnit);
  hm = fill3( ...
1.0005 * x(iM) + mRad * (bUnit(1) .* cos(iTh*pi/180) + cUnit(1) .* sin(iTh*pi/180)), ...
1.0005 * y(iM) + mRad * (bUnit(2) .* cos(iTh*pi/180) + cUnit(2) .* sin(iTh*pi/180)), ...
1.0005 * z(iM) + mRad * (bUnit(3) .* cos(iTh*pi/180) + cUnit(3) .* sin(iTh*pi/180)), ...
rColourx);
%    set(hm, 'EdgeColor', 0.95*( [1 1 1] - rCoef(iR) * ([1 1 1] - iColour) ) );
%    set(hm, 'FaceColor', 0.95*( [1 1 1] - rCoef(iR) * ([1 1 1] - rColouri) ) );
  hg = fill3( ...
x(iM) + gRad * (bUnit(1) .* cos(iTh*pi/180) + cUnit(1) .* sin(iTh*pi/180)), ...
y(iM) + gRad * (bUnit(2) .* cos(iTh*pi/180) + cUnit(2) .* sin(iTh*pi/180)), ...
z(iM) + gRad * (bUnit(3) .* cos(iTh*pi/180) + cUnit(3) .* sin(iTh*pi/180)), ...
rColourg);
end;


% Plot hemisphere joint
aUnit = [...
cos(hEl*pi/180) .* cos(hAz*pi/180); 
cos(hEl*pi/180) .* sin(hAz*pi/180);
sin(hEl*pi/180) ];
bUnit = [-sin(hAz*pi/180); cos(hAz*pi/180); 0];
cUnit = cross(aUnit, bUnit);
hh = fill3( ...
1.01 * r2 * (bUnit(1) .* cos(iTh*pi/180) + cUnit(1) .* sin(iTh*pi/180)), ...
1.01 * r2 * (bUnit(2) .* cos(iTh*pi/180) + cUnit(2) .* sin(iTh*pi/180)), ...
1.01 * r2 * (bUnit(3) .* cos(iTh*pi/180) + cUnit(3) .* sin(iTh*pi/180)), ...
aColour);


% Plot cable duct
xc1 = r2 * cos(cEl*pi/180) .* cos(cAz*pi/180);
yc1 = r2 * cos(cEl*pi/180) .* sin(cAz*pi/180);
zc1 = r2 * sin(cEl*pi/180);
aUnit = [...
cos(cEl*pi/180) .* cos(cAz*pi/180); 
cos(cEl*pi/180) .* sin(cAz*pi/180);
sin(cEl*pi/180) ];
bUnit = [-sin(cAz*pi/180); cos(cAz*pi/180); 0];
cUnit = cross(aUnit, bUnit);
hc1 = fill3( ...
1.0005 * xc1 + cRad1 * (bUnit(1) .* cos(iTh*pi/180) + cUnit(1) .* sin(iTh*pi/180)), ...
1.0005 * yc1 + cRad1 * (bUnit(2) .* cos(iTh*pi/180) + cUnit(2) .* sin(iTh*pi/180)), ...
1.0005 * zc1 + cRad1 * (bUnit(3) .* cos(iTh*pi/180) + cUnit(3) .* sin(iTh*pi/180)), ...
rColourx);
hc1 = fill3( ...
xc1 + cRad2 * (bUnit(1) .* cos(iTh*pi/180) + cUnit(1) .* sin(iTh*pi/180)), ...
yc1 + cRad2 * (bUnit(2) .* cos(iTh*pi/180) + cUnit(2) .* sin(iTh*pi/180)), ...
zc1 + cRad2 * (bUnit(3) .* cos(iTh*pi/180) + cUnit(3) .* sin(iTh*pi/180)), ...
rColourg);


% Plot fixings
for iF = 1:nF,
  xf = r2 * cos(fEl(iF)*pi/180) .* cos(fAz(iF)*pi/180);
  yf = r2 * cos(fEl(iF)*pi/180) .* sin(fAz(iF)*pi/180);
  zf = r2 * sin(fEl(iF)*pi/180);
  aUnit = [...
cos(fEl(iF)*pi/180) .* cos(fAz(iF)*pi/180); 
cos(fEl(iF)*pi/180) .* sin(fAz(iF)*pi/180);
sin(fEl(iF)*pi/180) ];
  bUnit = [-sin(fAz(iF)*pi/180); cos(fAz(iF)*pi/180); 0];
  cUnit = cross(aUnit, bUnit);
  hf = fill3( ...
1.0005 * xf + fRad1 * (bUnit(1) .* cos(iTh*pi/180) + cUnit(1) .* sin(iTh*pi/180)), ...
1.0005 * yf + fRad1 * (bUnit(2) .* cos(iTh*pi/180) + cUnit(2) .* sin(iTh*pi/180)), ...
1.0005 * zf + fRad1 * (bUnit(3) .* cos(iTh*pi/180) + cUnit(3) .* sin(iTh*pi/180)), ...
rColourx);
  hc1 = fill3( ...
xf + fRad2 * (bUnit(1) .* cos(iTh*pi/180) + cUnit(1) .* sin(iTh*pi/180)), ...
yf + fRad2 * (bUnit(2) .* cos(iTh*pi/180) + cUnit(2) .* sin(iTh*pi/180)), ...
zf + fRad2 * (bUnit(3) .* cos(iTh*pi/180) + cUnit(3) .* sin(iTh*pi/180)), ...
rColourf);
end;% for

axis(r2*[-1.2 1.2 -1.2 1.2 -1.2 1.2]);
drawnow;
hold off;
xlabel('x');
ylabel('y');
zlabel('z');

% orient portrait; print -dpdf coverage2.pdf