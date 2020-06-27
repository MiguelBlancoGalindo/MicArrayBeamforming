function [N_grid,array_pos] = elements_rectangular_grid_within_circle(N_squared_root,d,dx,dy,r)
% function [N_grid,array_pos] = elements_rectangular_grid_within_circle(N_squared_root,d,dx,dy,r)
% function that calculates the micrphone positions of a rectangular grid
% within a circle
%
% input parameters:
%   N_squared_root: number of microphones along a single dimension of an 
%       equivalent square array.
%   d: microphone spacing
%   dx: microphone spacing in the x dimension.
%   dy: microphone spacing in the x dimension.
%   r: radius of the circle where the array is to be confined. 
%
% output arguments: 
%   N_grid: number of microphones of the rectangular grid array. 
%   array_pos: microphone array positions in Cartesian coordinates. 

N_grid = 0;
indices = 0;
N_square = N_squared_root^2;

array_pos = zeros(N_grid,3);
array_pos1 = zeros(N_square,3);
distance = zeros(N_square,1);

for l=1:N_squared_root
    for m=1:N_squared_root
        array_pos1((l-1)*N_squared_root+m,1:2) = [-(l-1)*d (m-1)*d];
    end
end
%centering the array
array_pos1(:,1) = array_pos1(:,1) + max(abs(array_pos1(:,1)))/2 + dx;
array_pos1(:,2) = array_pos1(:,2) - max(abs(array_pos1(:,2)))/2 + dy;

for l=1:N_square
    distance(l) = sqrt(array_pos1(l,1)^2 + array_pos1(l,2)^2 + array_pos1(l,3)^2);
    if distance(l) <= r
        N_grid = N_grid + 1;
        indices(N_grid) = l;
    end
end

for n=1:N_grid
    array_pos(n,:) = array_pos1(indices(n),:);
end

end