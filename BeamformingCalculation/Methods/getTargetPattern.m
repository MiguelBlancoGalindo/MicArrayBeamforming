function d = getTargetPattern(Config)
% function d = getTargetPattern(Config)
% function that obtains the target directivity function from the settings
% specified in Config.
%
% input argument:
%   Config: configuration struct containing all the settings.
%
% output argument:
%   d: target directivity function. 

if isfield(Config.Filter,'ls') && isfield(Config.Filter.ls,'CurrMode') && isfield(Config.Filter.ls,'CurrOrder')
    method = 'ls';
elseif isfield(Config.Filter,'mn') && isfield(Config.Filter.mn,'CurrMode') && isfield(Config.Filter.mn,'CurrOrder')
    method = 'mn';
end
theta = Config.ArrayMan.LookAng.El;
phi = Config.ArrayMan.LookAng.Az;
L = max(length(theta),length(phi));
phi_t = Config.Filter.CurrTarget(2);
theta_t = Config.Filter.CurrTarget(1);
[~,~,itar] = getNearestAngle(Config,theta_t,phi_t);

if isfield(Config.Filter,'CurrInts') && ~isempty(Config.Filter.CurrInts)
intAngles = Config.Filter.CurrInts;
nints=length(intAngles)/2;
iint = zeros(nints,1);
    for kint=1:nints
        [~,~,iint(kint)] = getNearestAngle(Config,intAngles(kint*2-1),intAngles(kint*2));
    end
else
    iint=[];
end

if strcmpi(Config.Filter.(sprintf(method)).CurrMode,'hypercardioid')
    N=Config.Filter.(sprintf(method)).CurrOrder;
    if strcmp(Config.ArrayMan.SteerSpace,'3D') 
        [phi_P,theta_P,~] = mysph2physph(phi,theta,1);
        phi_P = pi/180.*phi_P;
        theta_P = pi/180.*theta_P;
        [phi_tP,theta_tP,~] = mysph2physph(phi_t,theta_t,1);
        phi_tP = pi/180.*repmat(phi_tP,L,1);
        theta_tP = pi/180.*repmat(theta_tP,L,1);
        cosTheta = cos(theta_tP).*cos(theta_P) + cos(phi_tP - phi_P) .* sin(theta_tP).*sin(theta_P);

    % Rafaely: Fundamentals of Spherical Arrays
        d=real(legendreP(N+1,cosTheta)-legendreP(N,cosTheta))./((N+1).*(cosTheta-1));
        d(itar) = 1;

    else
        if strcmp(Config.ArrayMan.SteerSpace,'2Del')
            omega = deg2rad(theta);
            omega_t = deg2rad(theta_t);
        elseif strcmp(Config.ArrayMan.SteerSpace,'2Daz')
            omega = deg2rad(phi);
            omega_t = deg2rad(phi_t);
        else
            error('currently not implemented in 3D');
        end
        F_ds=zeros(L,N+1);
        modeweights=[1,2*ones(1,N)];
        a=modeweights./sum(modeweights);

        % sum modal components of hypercardioids
        for n=0:N,F_ds(:,n+1)=a(n+1).*(cos(n*(omega-omega_t)));end
        d=sum(F_ds,2);

    end

elseif strcmp(Config.Filter.(sprintf(method)).CurrMode,'cardioid')
    N=Config.Filter.(sprintf(method)).CurrOrder;

    if strcmp(Config.ArrayMan.SteerSpace,'3D') 

        [phi_P,theta_P,~] = mysph2physph(phi,theta,1);
        phi_P = pi/180.*phi_P;
        theta_P = pi/180.*theta_P;
        [phi_tP,theta_tP,~] = mysph2physph(phi_t,theta_t,1);
        phi_tP = pi/180.*repmat(phi_tP,L,1);
        theta_tP = pi/180.*repmat(theta_tP,L,1);
        cosTheta = cos(theta_tP).*cos(theta_P) + cos(phi_tP - phi_P) .* sin(theta_tP).*sin(theta_P);

% Zotter and Frank (2019): "Ambisonics: A Practical 3D Audio Theory for
% Recording, Studio Production, Sound Reinforcement, and Virtual Reality"
        d=(1 + cosTheta).^N./2^N;

    else
        if strcmp(Config.ArrayMan.SteerSpace,'2Del')
            omega = deg2rad(theta);
            omega_t = deg2rad(theta_t);
        elseif strcmp(Config.ArrayMan.SteerSpace,'2Daz')
            omega = deg2rad(phi);
            omega_t = deg2rad(phi_t);
        end

        d=(1 + cos(omega-omega_t)).^N./2^N;
    end

elseif strcmpi(Config.Filter.(sprintf(method)).CurrMode,'chebyshev')
    % Nth order hypercardioid
    N = Config.Filter.(sprintf(method)).CurrOrder;
    BW = Config.Filter.(sprintf(method)).BW;
    BW = BW*pi/180;
    if strcmp(Config.ArrayMan.SteerSpace,'3D') 
    [phi_P,theta_P,~] = mysph2physph(phi,theta,1);
    phi_P = pi/180.*phi_P;
    theta_P = pi/180.*theta_P;
    [phi_tP,theta_tP,~] = mysph2physph(phi_t,theta_t,1);
    phi_tP = pi/180.*repmat(phi_tP,L,1);
    theta_tP = pi/180.*repmat(theta_tP,L,1);
    Theta = acos(cos(theta_tP).*cos(theta_P) + cos(phi_tP - phi_P) .* sin(theta_tP).*sin(theta_P));

    elseif strcmp(Config.ArrayMan.SteerSpace,'2Del')
        omega = deg2rad(theta);
        omega_t = deg2rad(theta_t);
    elseif strcmp(Config.ArrayMan.SteerSpace,'2Daz')
        omega = deg2rad(phi);
        omega_t = deg2rad(phi_t);
    end

    %chebyshev beampattern from "Spherical Microphone Array Beamforming"
    %Rafaely 2010
    alpha = omega-omega_t;
    alpha0 = BW/2;
    x0 = cos(pi/(2*N))./cos(alpha0/2);
    R = cosh(N*acosh(x0));
    x = x0*cos(alpha/2);
    d = 1/R.*chebyshevT(N,x);

elseif strcmp(Config.Filter.(sprintf(method)).CurrMode,'binary')
    % zeros apart from a one for the target angle
    d=zeros(size(ArrayMan(1).a,2),1);
    d(itar)=1;
elseif strcmp(Config.Filter.(sprintf(method)).CurrMode,'notch')
    error('notch mode not implemented');
elseif isempty(Config.Filter.(sprintf(method)).CurrMode)
    error(['no ' (sprintf(method)) '  target beampattern specified: check Config.Filter.' (sprintf(method)) '.CurrMode']);
else
    error(['no valid' (sprintf(method)) ' target beampattern: check Config.Filter.' (sprintf(method)) '.CurrMode']);
end

end

