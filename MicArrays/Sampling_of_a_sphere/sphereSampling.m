function SA = sphereSampling(SA,grid)

if strcmpi(grid,'Spk')
if ~isfield(SA,'SpkGrid')
    error(['inexistent field ''', 'SpkGrid' '''' ]);
else
    gridStruct = SA.SpkGrid;
    gridStr = 'SA.SpkGrid';
    number = 'L';
    
end

elseif strcmpi(grid,'Mic')
    if ~isfield(SA.Array,'MicGrid') || isempty(SA.Array.MicGrid)
        error(['inexistent field ''', 'MicGrid' '''' ]);
    else
        gridStruct = SA.Array.MicGrid;
        gridStr = 'SA.Array.MicGrid';
        number = 'M';
        
    end
else
    error(['incorret grid to sample the sphere on. Choose from ''', 'Spk', ''' or ''', 'Mic' '''' ]);
end

if ~isfield(gridStruct,'SamplingScheme') || isempty(gridStruct.SamplingScheme)
    error(['inexistent field ''', 'SamplingScheme', ''' in ''', gridStr '''' ]);
else
    sampling = gridStruct.SamplingScheme;
end
if ~isfield(gridStruct,'R') || isempty(gridStruct.R)
    error(['inexistent field ''', 'R', ''' in ''', gridStr '''' ]);
else
    r = gridStruct.R;
end

if strcmpi(sampling,'equiangular')
    if ~isfield(gridStruct,'N') || isempty(gridStruct.N)
        if ~isfield(gridStruct,number) || isempty(gridStruct.(sprintf(number)))
            if ~isfield(gridStruct,strcat(number,'phi')) || isempty(gridStruct.(sprintf(strcat(number,'phi'))))
                error(['inexistent field ''', 'N', ''' or ''', number, ''' in ''', gridStr '''' ]);
            else
                Lphi = gridStruct.(sprintf(strcat(number,'phi')));
                L = Lphi^2;
            end
        else
            L = gridStruct.(sprintf(number));
        end
        N = floor(sqrt(L/4)-1);
    else
        N = gridStruct.N;
    end
    theta = (1/2+(0:2*N+1)).*180/(2*N+2);
    phi = (0:2*N+1)*360/(2*N+2);
    alpha = zeros(2*N+2,1);
    for t=1:length(theta)
        temp = zeros(N+1,1);
        for q=0:N
            temp(q+1) = 1/(2*q+1)*sind((2*q+1)*theta(t));
        end
        alpha(t) = 2*pi/(N+1)^2*sind(theta(t))*sum(temp);
    end
    [theta,phi] = meshgrid(theta,phi);
    theta = reshape(theta,4*(N+1)^2,1);
    phi = reshape(phi,4*(N+1)^2,1);
    [alpha,~] = meshgrid(alpha,ones(2*N+2,1));
    alpha = reshape(alpha,4*(N+1)^2,1);
    [phi1,theta1,~] = physph2mysph(phi,theta,r);
    [cart(:,1),cart(:,2),cart(:,3)] = sph2car(phi1,theta1,r);
    L = 4*(N+1)^2;
elseif strcmpi(sampling,'Gaussian')
    if ~isfield(gridStruct,'N') || isempty(gridStruct.N)
        if ~isfield(gridStruct,number) || isempty(gridStruct.(sprintf(number)))
            if ~isfield(gridStruct,strcat(number,'phi')) || isempty(gridStruct.(sprintf(strcat(number,'phi'))))
                error(['inexistent field ''', 'N', ''' or ''', number, ''' in ''', gridStr '''' ]);
            else
                Lphi = gridStruct.(sprintf(strcat(number,'phi')));
                %N = Lphi/2-1;
                L = Lphi^2/2;
            end
        else
            L = gridStruct.(sprintf(number));
        end
        N = floor(sqrt(L/2)-1);
    else
        N = gridStruct.N;
    end
    theta = (1/2+(0:N)).*180/(N+1);
    phi = (0:2*N+1)*360/(2*(N+1));
    alpha = zeros(N+1,1);
    for t=1:length(theta)
        P = legendreP(N+2,cosd(theta(t)));
        alpha(t) = pi/(N+1)*2*(1-cosd(theta(t)).^2)/((N+2)^2.*P.^2);
    end
    [theta,phi] = meshgrid(theta,phi);
    theta = reshape(theta,2*(N+1)^2,1);
    phi = reshape(phi,2*(N+1)^2,1);
    [alpha,~] = meshgrid(alpha,ones(2*(N+1),1));
    alpha = reshape(alpha,2*(N+1)^2,1);
    [phi1,theta1,~] = physph2mysph(phi,theta,r);
    [cart(:,1),cart(:,2),cart(:,3)] = sph2car(phi1,theta1,r);
    L = 2*(N+1)^2;
elseif strcmpi(sampling,'uniform')
    
    if ~isfield(gridStruct,number) || isempty(gridStruct.(sprintf(number)))
        if ~isfield(gridStruct,'N') || isempty(gridStruct.N)
            error(['inexistent field ''', 'N', ''' or ''', number, ''' in ''', gridStr '''' ]);
        else
            N = gridStruct.N;
            L = (N+1)^2;
        end
    else
        L = gridStruct.(sprintf(number));
    end
    [cart,~] = ParticleSampleSphere('N',L);
    cart = cart*r;
    [phi,theta,~] = car2sph(cart(:,1),cart(:,2),cart(:,3));
    %converting to Physics coordinate system
    [phi,theta,~] = mysph2physph(phi,theta,r);
    phi(phi<0) = phi(phi<0)+360;
    alpha = repmat(4*pi/L,L,1);
    N = floor(sqrt(L)-1);
elseif strcmpi(sampling,'Lebedev')
    if ~isfield(gridStruct,number) || isempty(gridStruct.(sprintf(number)))
        if ~isfield(gridStruct,'N') || isempty(gridStruct.N)
            error(['inexistent field ''', 'N', ''' or ''', number, ''' in ''', gridStr '''' ]);
        else
            N = gridStruct.N;
            L = (N+1)^2;
        end
    else
        L = gridStruct.(sprintf(number));
    end
    
    [leb] = getLebedevSphere(L);
    cart = [leb.x leb.y leb.z]; 
    [phi,theta,~] = car2sph(leb.x,leb.y,leb.z);
    %converting to Physics coordinate system
    [phi,theta,~] = mysph2physph(phi,theta,r);
    phi(phi<0) = phi(phi<0)+360;
    alpha = leb.w;
    N = floor(sqrt(L)-1);
elseif strcmpi(sampling,'truncatedIcosahedron') 
    L = gridStruct.(sprintf(number));
    if isempty(L) 
        if isfield(gridStruct,'N') && gridStruct.N==4
            N = gridStruct.N;
            L = 32; 
            warning('32 sampling points for truncated Icosahedron to match up to 4th order');
        else
            error('number of sampling points or spherical harmonic order not specified');
        end
    elseif L~=32
        error(['Truncated Icosahedron not compatible with ' num2str(L) ' sampling points']);
    else
    load('EigenmikeManualSphCoor.mat');
    theta = EigenmikeManualSphCoor(:,1);
    phi = EigenmikeManualSphCoor(:,2);
    [phi1,theta1,~] = physph2mysph(phi,theta,r);
    [cart(:,1),cart(:,2),cart(:,3)] = sph2car(phi1,theta1,r);
    alpha = repmat(4*pi/L,L,1);
    N = floor(sqrt(L)-1);
    end
else
    error(['Unknown sampling scheme ' sampling]);
end

if strcmpi(grid,'Spk')
    SA.SpkGrid.SteerDir = [theta phi r.*ones(length(theta),1)];
    SA.SpkGrid.SteerDirCar = cart;
    SA.SpkGrid.QWeights = alpha;
    SA.SpkGrid.N = N;
    SA.SpkGrid.L = L;
elseif strcmpi(grid,'Mic')
    SA.Array.MicGrid.SteerDir = [theta phi r.*ones(length(theta),1)];
    SA.Array.MicGrid.SteerDirCar = cart;
    SA.Array.MicGrid.QWeights = alpha;
    SA.Array.MicGrid.N = N;
    SA.Array.MicGrid.M = L;
end
    
end
        