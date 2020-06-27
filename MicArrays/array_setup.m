function Config=array_setup(Config)
% function Config=array_setup(Config)
% function that calculates the microphone array positions based on the
% parameters specified in Config. 
%
% input and output arguments: 
%   Config: configuration structure containing all the settings.
%
%   Config.ArrayGeo: geometry of microphone array. Choose from: 
%       'linear_bs' (broadside), 'linear_ef' (endfire), 'circular', 
%       'dualcircular', 'multicircular', 'square', 'rectangular', 
%       'circular_grid', 'circular_rigid_sphere', 'circular_rigid_cylinder',
%       'stacked_circular' (vertically stacked circular), 
%       'stacked_circular_rigid_cylinder', 'spherical', 'spherical_rigid_sphere'

M = Config.Array.M; %number of microphones
r = Config.Array.r; %equivalent array radius
d = Config.Array.Spacing;  %microphone array spacing 
dim = Config.Array.Orientation; %dimension over which to span the array
NC = Config.Array.NCircles; %number of circles (if applicable)

Tri = 0;
MicPos = zeros(M,3); % this is a matrix with 3 columns corresponding to the 
% x, y, z coordinates of the microphones. 

MicPosSph = [];
MicPosCyl = [];
%%%%%%%%%%%%%%%%     linear array     %%%%%%%%%%%%%
if strcmp(Config.Array.Geo,'linear') || strcmp(Config.Array.Geo,'linear_bs') || strcmp(Config.Array.Geo,'linear_ef')
    if strcmp(Config.Array.Geo,'linear')  
        %line array is placed along dim axis with first sensor at x=d/2
        %sensors are numbered from right to left
        if strcmp(dim,'x')
            a=1;  
        elseif strcmp(dim,'y')
            a=2;
        elseif strcmp(dim,'z')
            a=3;
        else
            error('please introduce "x","y" or "z" as the dimension for the linear array to span over');
        end
    elseif strcmp(Config.Array.Geo,'linear_bs')
        a=1;  
    elseif strcmp(Config.Array.Geo,'linear_ef')
        a=2;
    end

    for n=1:M
        MicPos(n,a) = -(n-1)*d;
    end
    %centering the array
    MicPos(:,a) = MicPos(:,a) + max(abs(MicPos(:,a)))/2;

%%%%%%%%%%%%%%%%     circular array     %%%%%%%%%%%%%
elseif strcmp(Config.Array.Geo,'circular') || strcmp(Config.Array.Geo,'circular_rigid_sphere') || strcmp(Config.Array.Geo,'circular_rigid_cylinder') ||...
        strcmp(Config.Array.Geo,'circular rigid sphere') || strcmp(Config.Array.Geo,'circular rigid cylinder') 
    %sensors are numbered from the top-most sensor counterclockwise

    if isempty(r) && isempty(d)
        error('not enough input arguments');
    else
        if isempty(r)
            %the radius of the circumference is derived from the sensor
            %spacing and the angle beteween sensors (2*pi/N) using the
            %expression below
            %r=d/(sqrt(2)*sqrt(1-cos(2*pi/M)));
            r=d/(2*sin(pi/M));

        end

        %the vector a assings the sine and cosine values to the right
        %coordinate depending on the rotation
        
        %angle formed between the top-most sensor and the centre of the 
        % circle and the sensor of interest and the centre
        alpha = [360/M*(0:M-1)]'; 
        if strcmp(dim,'xy') || strcmp(dim,'yx')
            a=[1 2];  
            phi=alpha;             
            theta = zeros(M,1);
        elseif strcmp(dim,'yz') || strcmp(dim,'zy')
            a=[2 3];
            phi = [zeros(ceil(M/2),1); 180*ones(floor(M/2),1)];
            %restricting the values between 0 and 180
            theta = alpha -2.*floor(alpha./180).*mod(alpha,180);
            %converting from zenith to elevation
            theta = 90 - theta;
        elseif strcmp(dim,'zx') || strcmp(dim,'xz')
            a=[1 3];
            phi = [90*ones(ceil(M/2),1); -90*ones(floor(M/2),1)];
            %restricting the values between 0 and 180
            theta = alpha -2.*floor(alpha./180).*mod(alpha,180);
            %converting from zenith to elevation
            theta = 90 - theta;
        else
            error(['please introduce "xy","yz" or "xz" as the dimensions'...
            'for the linear array to span over']);
        end
        alpha = alpha*pi/180;
        MicPosSph = [phi,theta,repmat(r,M,1)];
        MicPos(:,a) = [r*sin(alpha) r*cos(alpha)];

    end
%%%%%%%%%%%%%%%%     circular array     %%%%%%%%%%%%%
elseif strcmp(Config.Array.Geo,'stacked_circular') || strcmp(Config.Array.Geo,'stacked_circular_rigid_cylinder')
    %sensors are numbered from the top-most sensor counterclockwise

    if isempty(r) && isempty(d)
        error('not enough input arguments');
    else
        if isempty(r)
            %the radius of the circumference is derived from the sensor
            %spacing and the angle beteween sensors (2*pi/N) using the
            %expression below
            %r=d/(sqrt(2)*sqrt(1-cos(2*pi/M)));
            r=d/(2*sin(pi/M));

        end
        if isempty(NC)
            %calculating the maximum number of circles that is a divisor of
            %M and is smaller than the resulting number of mikes per circle
            NC=find(M/divisors(M)>M/divisors(M),'last');
        else
            if mod(M,NC)~=0
                error('The number of microphones is not divisible by the number of circles')
            end
        end
        Mint = floor(M/NC);
        Mleft = M-Mint*NC;
        Mc = repmat(Mint,NC,1);
        Mc(1:Mleft) = Mc(1:Mleft) + 1;
        im = 1;
        
        if strcmp(Config.Array.MultiCircleAlignment,'notaligned')
            phi = zeros(M,1);
            for ic = 1:NC
                %angle formed between the top-most sensor and the centre of the 
                % circle and the sensor of interest and the centre
                %each circle is offset by and even amount
                phi(im:Mc(ic)+im-1) = (360/Mc(ic)*(0:Mc(ic)-1)).' + 360/Mint/NC*(ic-1); 
                im = im + Mc(ic);
            end
        else
            phi = repmat([360/Mint*(0:Mint-1)]',NC,1);
        end
        z = zeros(M,1);
        im = 1 + Mc(1);
        for ic=2:NC
            z(im:Mc(ic)+im-1) = (ic-1)*Config.Array.StackedCircleSpacing;
            im = im + Mc(ic);
        end
        %centering array vertically
        z = z - max(abs(z))/2;
         
        if strcmp(dim,'xy') || strcmp(dim,'yx')
            
            MicPosCyl = [phi,z,repmat(r,M,1)];
            [MicPos(:,1),MicPos(:,2),MicPos(:,3)] = cyl2car(phi,z,repmat(r,M,1));

        else
            error(['please introduce "xy","yz" or "xz" as the dimensions'...
            'for the stacked circular array to span over']);
        end

    end
%%%%%%%%%%%%%%%%     dualcircular array     %%%%%%%%%%%%%
elseif strcmp(Config.Array.Geo,'dualcircular')  
    %sensors are numbered from the right-most sensor counterclockwise
    %spacing d corresponds to the microphone spacing of the
    %outermost circular array.

    if isempty(r) && isempty(d)
        error('not enough input arguments');
    else
        %If the number of sensors is even each circumference has the
        %same number. Otherwise the outermost circumference has one
        %sensor more (N1)
        if mod(M,2)>0
            M1=ceil(M/2);
        else
            M1=M/2;
        end
        M2=M-M1;
        if isempty(r)          
            
            %the radius of the outermost circumference is derived from the sensor
            %spacing and the angle beteween sensors (2*pi/N1) using the
            %expression below
            r=zeros(1,2);
            r(2)=d/(2*sin(pi/(M/2)));
            %r(2)=d/(sqrt(2)*sqrt(1-cos(2*pi/N1)));
            %radius of the innermost circumference 
                        
            %calculating the inner circle with the same ratio as the Phil's
            %dual-circular array, and assuming the spacing is kept for the
            %inner
            r(1) = r(2)*0.083/0.107;
            
%             %The radial distance between both circumferences is given by 
%             delta_r=input(['Introduce the radial distance (in meters) '...
%                 'between both circular arrays ']);
%             r(1)=r(2)-delta_r;
        end
        if length(r)~=2
            error(['the radius value must be a 1x2 vector with the two'...
                ' microphone radii']);
        end
        if strcmp(Config.Array.MultiCircleAlignment,'notaligned')
            alpha_offset = 360/M1/2;
        else
            alpha_offset = 0;
        end
        alpha1 = [360/M1*(0:M1-1)]';
        alpha2 = [360/M2*(0:M2-1)]' + alpha_offset;
        alpha = [alpha1; alpha2];
        %the vector a assings the sine and cosine values to the right
        %coordinate depending on the rotation

        if strcmp(dim,'xy') || strcmp(dim,'yx')
            a=[1 2];  
            phi=alpha;             
            theta = zeros(M,1);
        elseif strcmp(dim,'yz') || strcmp(dim,'zy')
            a=[2 3];
            phi = [zeros(ceil(M1/2),1); 180*ones(floor(M1/2),1);...
                zeros(ceil(M2/2),1); 180*ones(floor(M2/2),1)];
            %restricting the values between 0 and 180
            theta = alpha -2.*floor(alpha./180).*mod(alpha,180);
            %converting from zenith to elevation
            theta = 90 - theta;
        elseif strcmp(dim,'zx') || strcmp(dim,'xz')
            a=[1 3];
            phi = [90*ones(ceil(M/2),1); -90*ones(floor(M/2),1)];
            %restricting the values between 0 and 180
            theta = alpha -2.*floor(alpha./180).*mod(alpha,180);
            %converting from zenith to elevation
            theta = 90 - theta;
        else
            error(['please introduce "xy","yz" or "xz" as the dimensions'...
            'for the linear array to span over']);
        end
        alpha1 = alpha1*pi/180;
        alpha2 = alpha2*pi/180;
        MicPosSph = [phi,theta,[repmat(r(2),M1,1);repmat(r(1),M2,1)]];
        MicPos(:,a) = [r(2)*sin(alpha1) r(2)*cos(alpha1); r(1)*sin(alpha2) r(1)*cos(alpha2)];
        
    end
    
%%%%%%%%%%%%%%%%     multicircular array     %%%%%%%%%%%%%    
elseif strcmp(Config.Array.Geo,'multicircular') 
    theta = 0;
    phi = zeros(M,1);
    rn = r*(NC-([1:NC]'-1))/NC;    %radius of each circumferece
    r = zeros(M,1);     %reseting r to be the value to be saved
    if strcmp(Config.Array.MultiCircleAlignment,'notaligned')
        d0 = 2*pi*sum(rn)/M;
        Mn = round(2*pi*rn./d0);
        iM = 1;
        for iC = 1:NC
            alpha = 360/Mn(iC)*[1:Mn(iC)]';
            phi(iM:iM+Mn(iC)-1) = alpha;
            r(iM:iM+Mn(iC)-1) = repmat(rn(iC),Mn(iC),1);
            iM = iM + Mn(iC);
        end
    else
        Mn = floor(M/NC);
        phi = repmat(360/Mn*[1:Mn]',NC,1);

        for iC = 1:NC
            r((iC-1)*Mn+1:iC*Mn) = repmat(rn(iC),Mn,1);
        end
        if mod(M,NC)==1 
            phi = [phi; 0];
            r = [r; 0];
        elseif mod(M,NC)~=1 && mod(M,NC)~=0
            error('number of mikes needs to be multiple of number of circumferences for "aligned" sensors');
        end
    end
    MicPosSph = [phi repmat(theta,M,1) r];
    [MicPos(:,1),MicPos(:,2),MicPos(:,3)] = sph2car(phi,repmat(theta,M,1),r);
   
            
%%%%%%%%%%%%%%%%     square array     %%%%%%%%%%%%%
elseif strcmp(Config.Array.Geo,'square')           
    if(sqrt(M)==floor(sqrt(M))) %N must be a perfect square
        %vector a assings the sine and cosine values to the right
        %coordinate depending on the rotation
        if strcmp(dim,'xy') || strcmp(dim,'yx')
            a=[1 2];  
        elseif strcmp(dim,'yz') || strcmp(dim,'zy')
            a=[3 2];
        elseif strcmp(dim,'zx') || strcmp(dim,'xz')
            a=[3 1];
        else
            error(['please introduce "xy","yz" or "xz" as the dimensions'...
            'for the planar array to span over']);
        end

        for l=1:sqrt(M)
            for m=1:sqrt(M)
                MicPos((l-1)*sqrt(M)+m,a) = [-(l-1)*d (m-1)*d];
            end
        end
    else
        error(['number of sensors must be a perfect square number'...
            'i.e. whose square root is exact']);
    end
    %centering the array
    MicPos(:,1) = MicPos(:,1) + max(abs(MicPos(:,1)))/2;
    MicPos(:,2) = MicPos(:,2) - max(abs(MicPos(:,2)))/2;
    
%%%%%%%%%%%%%%%%     rectangular array     %%%%%%%%%%%%%
elseif strcmp(Config.Array.Geo,'rectangular')
    if sqrt(M)~=floor(sqrt(M))
        My=1;
        Mx=M;
        exact = zeros(M,1);
        while My<Mx
            Mx=floor(M/My);
            if mod(M,My)==0
                exact(My)=My;
            end
            My=My+1;
        end
        My=max(exact);
        Mx=M/My;
    else
        My=sqrt(M);
        Mx=My;
    end

    if strcmp(dim,'xy') || strcmp(dim,'yx')
        a=[1 2];  
    elseif strcmp(dim,'yz') || strcmp(dim,'zy')
        a=[3 2];
    elseif strcmp(dim,'zx') || strcmp(dim,'xz')
        a=[3 1];
    else
        error(['please introduce "xy","yz" or "xz" as the dimensions'...
        'for the planar array to span over']);
    end

    for l=1:Mx
        for m=1:My
            MicPos((l-1)*My+m,a) = [-(l-1)*d (m-1)*d];
        end
    end

    %centering the array
    MicPos(:,1) = MicPos(:,1) + max(abs(MicPos(:,1)))/2;
    MicPos(:,2) = MicPos(:,2) - max(abs(MicPos(:,2)))/2;

%%%%%%%%%%%%%%%%     circular grid array     %%%%%%%%%%%%%
elseif strcmp(Config.Array.Geo,'circular_grid') || strcmp(Config.Array.Geo,'circular grid')
    
    [MicPos,circle_pos] = generate_rectangular_grid_within_circle(M, r, dim);
    
    figure;
    for n=1:size(MicPos,1)
        scatter(MicPos(n,2),MicPos(n,1),'ok'); hold on
    end
    for n=1:length(circle_pos)
        scatter(circle_pos(n,2),circle_pos(n,1),'.r'); hold on
    end
    axis equal
        
%%%%%%%%%%%%%%%%     spherical array     %%%%%%%%%%%%%
elseif strcmp(Config.Array.Geo,'spherical') || strcmp(Config.Array.Geo,'spherical_rigid_sphere') || strcmp(Config.Array.Geo,'spherical rigid sphere')
    %a nearly uniform distribution of the sensors on a sphere is 
    %carried out by this external function, by an optimsisation process
    %depending on the number of sensors. 
    %the order of the sensors around the sphere is unknown
    % Tri is the M-by-3 list of face-vertex connectivities
    if M==32
        load('EigenmikeManualSphCoor.mat');
        [phi,theta,r] = physph2mysph(EigenmikeManualSphCoor(:,2),EigenmikeManualSphCoor(:,1),r);
        [x,y,z]=sph2car(phi,theta,r);
        MicPos = [x,y,z];
        MicPosSph = [phi,theta,repmat(r,M,1)];
    elseif M==8
        MicPos = r/sqrt(3)*[-1 -1 -1; -1 -1 1; -1 1 -1; -1 1 1; ...
                    1 -1 -1; 1 -1 1; 1 1 -1; 1 1 1];
        [MicPosSph(:,1),MicPosSph(:,2),MicPosSph(:,3)] = car2sph(MicPos(:,1),MicPos(:,2),MicPos(:,3));
    elseif M==4
        MicPos = r/sqrt(3)*[-1 -1 -1; -1 1 1; 1 -1 1; 1 1 -1];
        [MicPosSph(:,1),MicPosSph(:,2),MicPosSph(:,3)] = car2sph(MicPos(:,1),MicPos(:,2),MicPos(:,3));
    else
        [MicPos,Tri] = ParticleSampleSphere('N',M);
        MicPos = MicPos*r;
        [MicPosSph(:,1),MicPosSph(:,2),MicPosSph(:,3)] = car2sph(MicPos(:,1),MicPos(:,2),MicPos(:,3));
    end


%         %other method
%         %[x,y,z,avgr] = psphere(n,demo);
%         demo=1; %displays the point on the sphere for each iteration
%         %avgr is the minimum distance between two points
%         [x,y,z,avgr] = psphere(N,demo);
%         MicPos(:,1) = x;
%         MicPos(:,2) = y;
%         MicPos(:,3) = z;
%         


%     figure
%     scatter(MicPos(:,1),MicPos(:,2));
%     axis 'equal';

elseif strcmp(Config.Array.Geo,'spiral') 
    [MicPos,MicPosCyl] = makeSpiralArray(M,r,Config.Array.L);
    
elseif strcmp(Config.Array.Geo,'ref_omni') 
    MicPos = zeros(1,3);
else
    error(['The array geometry is incorrect. Choose from ''', 'linear' ''', ''',...
    'circular' ''', ''', 'dualcircular' ''', ''', 'multicircular' ''', ''', ...
    'rectangular' ''', ''', 'square' ''', ''', 'spherical' ''', ''',...
    'circular_rigid_cylinder' ''', ''', 'circular_rigid_sphere' ''', ''',...
    'spherical_rigid_sphere' ''', ''','circular_grid',''' or ''','ref_omni', '''']);
end

Config.Array.MicPos = MicPos;
Config.Array.MicPosSph = MicPosSph;
Config.Array.MicPosCyl = MicPosCyl;
end


