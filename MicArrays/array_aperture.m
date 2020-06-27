function aperture = array_aperture(array_pos,array_geo, method)
% function aperture = array_aperture(array_pos,array_geo, method)
% function that calculates the (maximum) distance between all microphones.
%
% input arguments:
%   array_pos: microphone array positions in Cartesian coordinates (M x 3) 
%       for M microphones. 
%   array_geo: microphone array geometry.
%   method: type of calculation for rigid scatterers. Choose from
%   'near_field' or 'far_field'.
% 
% output arguments: 
%   aperture: maximum aperture of the array.

nMics = size(array_pos,1);
% the total number of distances to be calculated is equal to the number of
% combinations from nMics taken in pairs without repetition.
% combinations = factorial(nMics)/(factorial(2)*factorial(nMics-2));
%or equivalently for a more efficient implementation
combinations = nMics*(nMics-1);

aperture = zeros(combinations, 1);
count=1;

for n=1:nMics
    for m=n+1:nMics
        aperture(count) = sqrt((array_pos(n,1)-array_pos(m,1)).^2 + ...
            (array_pos(n,2)-array_pos(m,2)).^2 + ...
            (array_pos(n,3)-array_pos(m,3)).^2);
        count = count +1;
    end
end

aperture = max(aperture);

if strcmpi(array_geo,'spherical_rigid_sphere') || strcmp(array_geo,'circular_rigid_cylinder') || ...
    strcmpi(array_geo,'circular_rigid_sphere') || strcmp(array_geo,'stacked_circular_rigid_cylinder')
   
    centre = [mean(array_pos(:,1)) mean(array_pos(:,2)) mean(array_pos(:,3))];
    array_pos = array_pos - repmat(centre,nMics,1);
    [az,el,r] = car2sph(array_pos(:,1),array_pos(:,2),array_pos(:,3));
    r=mean(r);

    % calculating the distance as the sqrt of the sum of the squared 
    % distances of the arc on the horizontal plane and the z coordinate.

    d = zeros(nMics,nMics);

    if strcmpi(array_geo,'spherical_rigid_sphere')
        for m=1:nMics-1
            for m_prime = m+1:nMics
                d(m,m_prime) = real(r.*acos(array_pos(m,:)*array_pos(m_prime,:).'./r.^2));
            end
        end
    else
        for m=1:nMics-1
            for m_prime = m+1:nMics
                d(m,m_prime) = sqrt((r.*acos(array_pos(m,1:2)*array_pos(m_prime,1:2).'./r.^2)).^2 + (array_pos(m,3)-array_pos(m_prime,3)).^2);
            end
        end
    
    end

    aperture = max(max(d));
            
    %if microphones are uniformly distributed the distance will always be
    %d = pi*r;
        
    if strcmpi(method,'far_field')

        %since all baffled arrays are circular or spherical, the effective
        %aperture can be calculated from any microphone, e.g. the first one
        phi = az(1);
        theta = el(1);

        %free field travel path difference with respect to the centre for the look
        %direction phi=0, theta=0 (this one is used as a microphone will
        %presumably match this direction)
        free_field_path_diff = (array_pos(:,1)*sin(phi*pi/180) + array_pos(:,2)*cos(phi*pi/180))*cos(theta*pi/180) + array_pos(:,3)*sin(theta*pi/180); 
        path_diff = free_field_path_diff;
        for n=1:nMics
            if free_field_path_diff(n) < 0 
                path_diff(n) = -(d(1,n) - r*pi/2);
            end
        end
        
        aperture = max(path_diff)-min(path_diff);

    end

end
