% rs is the vector associated with the position of the sound source with
% respect to the origin of coordinates in spherical coordinates (azimuth
% and elevation only for plane waves)

% r is the Mx3 matrix associated with the positions of the microphones on
% the sphere in cartesian coordinates


function p = pressure_on_rigid_sphere_from_plane_wave(r,rs_hat,k)

M = size(r,1);
r_hat_norm2 = zeros(M,1);
for imics=1:M
    r_hat_norm2(imics) = sqrt(r(imics,1)^2 + r(imics,2)^2 + r(imics,3)^2);
end
r_hat_norm2 = median(r_hat_norm2);
r_hat = r/r_hat_norm2;

x=k*r_hat_norm2;
N = ceil(1.1*max(k)*r_hat_norm2);

%creating a matrix of the same size as r_hat
rs_hat = repmat(rs_hat,M,1);
Phi = dot(r_hat,rs_hat,2);  %angle difference between mic and source positions

p = 0;
%ptot = zeros(M,length(k),N+1);
for n=0:N
    %disp(['Spherical Harmonic Order: ' num2str(n)]);
    %spherical bessel function of first kind
    %jn = sqrt(pi./(2*x)).*besselj(n+1/2,x);
    
    %Spherical Hankel function of second kind
    %hn = sqrt(pi./(2*x)).*(besselj(n+1/2,x) + 1i*bessely(n+1/2,x));
    hn = sqrt(pi./(2*x)).*besselh(n+1/2,2,x);
    hn_1 = sqrt(pi./(2*x)).*besselh(n-1/2,2,x);
    %First derivative of Bessel function of first kind
    %djn = sqrt(pi./(2*x)).*besselj(n-1/2,x) - (n+1)./x.*jn;

    %First derivative of Hankel function of second kind
    %dhn = sqrt(pi./(2*x)).*(besselj(n-1/2,x) + 1i*bessely(n-1/2,x)) - (n+1)./x.*hn;
    dhn = hn_1 - (n+1)./x.*hn;
    
    %mode strength. It is the element that varies between open and rigid spheres
    %bn = jn-djn./dhn.*hn;
    %Wronskian relation
    %bn = 1i./(dhn.*x.^2); 
    bn = 1./dhn; 
    
    %Legendre polynomial
    Pn = legendreP(n,Phi);
    p = p + (1i).^n*Pn*(bn.*(2*n+1));
    %ptot(:,:,n+1) = (1i).^n*Pn*(bn.*(2*n+1));

end

%p = 1./(1i.*x.^2).*squeeze(sum(ptot,3));
p = p./(1i.*x.^2);

end
