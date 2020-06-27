function [d,dist]=calculate_array_spacing(varargin)
% function [d,dist]=calculate_array_spacing(varargin)
% function that calculates the (minimum) sensor spacing of a microphone array
%
% input arguments: choose from the following strings followed by a coma 
% and its value
%   'ArrayPos': microphone array positions in Cartesian coordinates. 
%   'ArrayGeo': microphone array geometry.
%   'mics': number of microphones.
%   'r': radius of the equivalent maximum aperture '2r'. 
%
% output arguments: 
%   d: minimum microphone spacing.
%   dist: all microphone distances of adjacent pairs. 

minInputs = 1;
maxInputs = Inf;
narginchk(minInputs, maxInputs);

array_pos_n = 0;
array_geo_n = 0;
M_n = 0;
r_n = 0;
dist = [];
for n=1:nargin
    if strcmp(varargin{n},'ArrayPos')
        array_pos = varargin{n+1};
        array_pos_n = n+1;
        M = size(array_pos,1);
        break;
    end
end

for n=1:nargin
    if strcmp(varargin{n},'ArrayGeo')
        array_geo = varargin{n+1};
        array_geo_n = n+1;
        break;
    end
end

for n=1:nargin
    if strcmp(varargin{n},'mics')
        M = varargin{n+1};
        M_n = n+1;
        break;
    end
end

for n=1:nargin
    if strcmp(varargin{n},'r')
        r = varargin{n+1};
        r_n = n+1;
        break;
    end
end

if array_geo_n ~= 0 %if array_geo has been passed as a parameter
    if strcmp(array_geo,'circular_rigid_sphere') || strcmp(array_geo,'circular_rigid_cylinder') || ...
	strcmp(array_geo,'stacked_circular_rigid_cylinder')
        if array_pos_n==0 && (M_n==0 || r_n==0)
            error(['specify "mics" and "r" or "array_pos" as input parameters for the ' array_geo ' array']);
        elseif array_pos_n~=0
            centre = [mean(array_pos(:,1)) mean(array_pos(:,2)) mean(array_pos(:,3))];
            array_pos = array_pos - repmat(centre,M,1);
            r = zeros(M,1);
            %addpath '../Evaluation/sphere_evaluation_error/sphere_Miguel'
            for n=1:M 
                r(n) = sqrt(array_pos(n,1)^2 + array_pos(n,2)^2 + array_pos(n,3)^2);
            end
            r=mean(r);
        end
        if strcmp(array_geo,'stacked_circular_rigid_cylinder')
            NC = 1;
            z0 = array_pos(1,3);
            for n=1:M
                z=array_pos(n,3);
                if z~=z0
                    NC = NC +1;
                    z0 = z;
                end
            end
            d = 2*pi/(M/NC)*r;
        else
            d = 2*pi/M*r;
        end

    elseif strcmp(array_geo,'spherical') || strcmp(array_geo,'spherical_rigid_sphere') 
        if strcmp(array_geo,'spherical')
            sphere_type = 'open';
        else
            sphere_type = 'solid';
        end
        if array_pos_n==0 && (M_n==0 || r_n==0)
            error(['specify "mics" and "r" or "array_pos" as input parameters for the ' array_geo ' array']);
        elseif array_pos_n==0 && (M_n~=0 || r_n~=0)
            array_pos = array_setup_v7(array_geo, M, [], r);
        else
            centre = [mean(array_pos(:,1)) mean(array_pos(:,2)) mean(array_pos(:,3))];
            array_pos = array_pos - repmat(centre,M,1);
            r = zeros(M,1);
            %addpath '../Evaluation/sphere_evaluation_error/sphere_Miguel'
            for n=1:M 
                r(n) = sqrt(array_pos(n,1)^2 + array_pos(n,2)^2 + array_pos(n,3)^2);
            end
            r=mean(r);
        end
            [d_closest_mic,d_several_mics,~,allDist]=array_pos_deviation(array_pos,r,sphere_type);
            % spacing is taken as the average of all the distances
            % with respect to closest microphone.
            d = d_closest_mic;

    else
        if array_pos_n == 0 %if array_pos has not been passed as a parameter
            array_pos = array_setup_v7(array_geo,M,[],r);
        end
        %spacing calculated as min distance between consecutive elements
        for m=1:M-1
            dist(m) = sqrt((array_pos(m,1)-array_pos(m+1,1))^2 + ...
                (array_pos(m,2)-array_pos(m+1,2))^2 + ...
                (array_pos(m,3)-array_pos(m+1,3))^2);           
        end
        d = min(dist);
    end
end

        
