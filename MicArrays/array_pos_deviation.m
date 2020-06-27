function [avg_closest,avg,sigma,min_dis_vector]=array_pos_deviation(array_coor,r,sphere_type)
% function [avg_closest,avg,sigma,min_dis_vector]=array_pos_deviation(array_coor,r,sphere_type)
% This function calculates the mean and standard deviation of the relative
% distance of the different sampling points on the surface of a sphere
%
% input parameters:
%   array_coor: N by 3 matrix, whose columns are the x,y,z coordinates
%       of all the sensors of the sphere.
%   r: sphere radius.
%   sphere_type: string of characters which if 'solid' calculates the
%       distance of a solid sphere using the arc length between two sensors, 
%       and if 'open' calculates the euclidian distance between sensors 
%       corresponding to the actual distance for a wave to travel between
%       points.
%
% output arguments: 
%   avg_closest: value of the closest microphone distance averaged over all
%       sensors.
%   avg: average of all the 5/6 closest microphones averaged over all
%       sensors.
%   sigma: standard deviation of the distances of the 5/6 closest microphones. 
%   min_dis_vector: distances of the 5/6 closest microphones. 

N=size(array_coor,1);
euc_distance = 1000*ones(N-1,N-1);

for i=1:N-1
    for j=i+1:N
        euc_distance(i,j) = sqrt((array_coor(i,1)-array_coor(j,1)).^2 + ...
            (array_coor(i,2)-array_coor(j,2)).^2 + ...
            (array_coor(i,3)-array_coor(j,3)).^2);
    end
end

if strcmp(sphere_type,'solid')
    arc_distance = 2*r*asin(euc_distance./(2*r));
    distance = arc_distance;
elseif strcmp(sphere_type,'open')
    distance = euc_distance;
end

min_dis = nan(N-1,6);
%min_dis_pos = zeros(N-1,6);
distance_mod = distance;
count = zeros(N,1);
for i=1:N-1
    
    while count(i)<6
        
        [min_dist, min_dis_pos] = min(distance_mod(i,:));
        if min_dist~=1000
            count(i)=count(i)+1;
            min_dis(i,count(i)) = min_dist;
            count(min_dis_pos) = count(min_dis_pos)+1;
            distance_mod(i,min_dis_pos) = 1000;
            %if the 6th closest point is larger than 1.1 times the distance
            %to the 5th closest point then the 6th point will not be considered
            %as it is assumed that that point will only have 5 connections
            if count(i)==6 && min_dist > 1.1*min_dis(i,5)
                min_dis(i,count(i)) = nan;
            end
        else
            count(i) = 6;
        end
           
    end
end


count=0;
for i=1:N-1
    for j=1:6
        if ~isnan(min_dis(i,j))
            count = count+1;
            min_dis_vector(count) = min_dis(i,j);
%             if min_dis(i,j)< closest_mic_dis(i)
%                 closest_mic_dis(count) = min_dis(i,j);
%             end
        end
    end
end

%value of the closest microphone distance averaged over all sensors
avg_closest = mean(min(min_dis,[],2,'omitnan'),'omitnan');    
%average of all the 5/6 closest microphones averaged over all sensors
avg = mean(min_dis_vector);
sigma = std(min_dis_vector,1);

end
  
      
        
