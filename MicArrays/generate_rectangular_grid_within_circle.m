function [rect_grid_pos,circle_pos] = generate_rectangular_grid_within_circle(N, r, dim)
N_squared_root = 12;
N_square = N_squared_root^2;

d = 0.1;    %fixed grid which will be resized after
r1=d+d*0.5;  %initial radius to give 1 microphone in the area of the circle
%vector a assings the sine and cosine values to the right
%coordinate depending on the rotation
if nargin<5 || strcmp(dim,'xy') || strcmp(dim,'yx')
    a=[1 2];  
elseif strcmp(dim,'yz') || strcmp(dim,'zy')
    a=[3 2];
elseif strcmp(dim,'zx') || strcmp(dim,'xz')
    a=[3 1];
else
    error(['please introduce "xy","yz" or "xz" as the dimensions'...
    'for the planar array to span over']);
end

%calculating circle within which the rectangular grid has to lie
N1=100;
circle_pos = zeros(N1,3);
for n=1:N1
    alpha=2*pi/N1*(n-1);  
    circle_pos(n,a) = [r*sin(alpha) r*cos(alpha)];
end

found=0;
largeN = 0;
dx = 0;
dy = 0;
fract = 1/20;

[N_grid,~] = elements_rectangular_grid_within_circle(N_squared_root,d,dx,dy,r1);
while N_grid<N 
    r1=r1+d*fract;
    [N_grid,grid_pos] = elements_rectangular_grid_within_circle(N_squared_root,d,dx,dy,r1);
end

if N_grid==N
    found = 1;
end

while found==0
    d = 0.1;    %fixed grid which will be resized after
    r1=d+d*0.5;  %initial radius to give 1 microphone in the area of the circle
    dx=d/2;
    [N_grid,~] = elements_rectangular_grid_within_circle(N_squared_root,d,dx,dy,r1);
    while N_grid<N 
        r1=r1+d*fract;
        [N_grid,grid_pos] = elements_rectangular_grid_within_circle(N_squared_root,d,dx,dy,r1);
    end
    if N==N_grid
        found=1;
    end
    if found == 0 && dy==0
        dy=d/2;
    elseif found == 0 && dy==d/2
        disp(['grid not found for N=' num2str(N)]);
        break;
    end
        
end

%resizing grid
if found==1
    ratio = r/r1;
    rect_grid_pos = grid_pos*ratio;
    %checking again with respect to the aperture
    aperture_rect_grid = array_aperture_v2(rect_grid_pos,'circular_grid');
    rect_grid_pos = rect_grid_pos*r/(aperture_rect_grid/2);
end

end