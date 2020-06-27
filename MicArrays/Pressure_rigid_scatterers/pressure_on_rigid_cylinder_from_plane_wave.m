% rs is the vector associated with the position of the sound source with
% respect to the origin of coordinates

% r is the vector associated with the position of the centre of the array 
% with respect to the origin of coordinates


function p = pressure_on_rigid_cylinder_from_plane_wave(r,rs_hat,k)

M = size(r,1);

[phi,theta,r_hat_norm2] = car2sph(r(:,1),r(:,2),r(:,3));
z=r(:,3);
%phi = phi*pi/180;
r_hat_norm2 = median(r_hat_norm2);
r_hat = r/r_hat_norm2;

[phis,thetas,rs_hat_norm2] = car2sph(rs_hat(:,1),rs_hat(:,2),rs_hat(:,3));

%angles are transformed to the standard physics coordinate system
[phi_P,theta_P,~] = mysph2physph(phi,theta,r_hat_norm2);
[phis_P,thetas_P,~] = mysph2physph(phis,thetas,r_hat_norm2);

phi_P = phi_P*pi/180;
theta_P = theta_P*pi/180;
phis_P = phis_P*pi/180;
thetas_P = thetas_P*pi/180;

x=k*r_hat_norm2;
N = ceil(1.1*max(k)*r_hat_norm2);
krsinthetas=k.*r_hat_norm2.*sin(thetas_P);

%phis_P = repmat(phis_P,M,1);
thetas_P = repmat(thetas_P,M,1);

e1 = exp(1i.*cos(thetas_P).*z*k);

%for the derivative of the Hankel function
% epsilon = [2*ones(1,N) 1 2*ones(1,N)];
% p = 0;
% p=zeros(M,length(k));
% pinc = zeros(N,1,M,length(k));
% pscat = zeros(N,N,M,length(k));
ptot = zeros(M,length(k),2*N+1);
%for imic=1:M
for n=-N:1:N
    in=n+N+1;
    %disp(['Spherical Harmonic Order: ' num2str(n)]);
    
    %Bessel function of first kind
    %Jn = besselj(n,x);

    %Hankel function of first kind
    %Hn = besselj(n+1/2,x) + 1i*bessely(n+1/2,x));
    %Hn = besselh(n,1,x);
    %H2nminus1 = besselh(n-1,2,x*cos(thetas));
    %H2nplus1 = besselh(n+1,2,x*cos(thetas));

    %First derivative of Bessel function of first kind
    %dJn = 1/2.*(besselj(n-1,x) - besselj(n+1,x));

    %First derivative of Hankel function of first kind
    %dHn = 1/2.*(besselh(n-1,x) - besselh(n+1,x));
    %First derivative of Hankel function of second kind
    %dH2n = 2*(H2nminus1-H2nplus1);
    %approximation for small values of x (not accurate around 0 for
    %numerical error I guess and as kr approaches 1)
%     factorialn=(sign(n)).^n.*factorial(abs(n));
%     dHn = 1i*factorialn./(pi*epsilon(n+N+1)).*(repmat(2,1,length(k))./x).^(n+1);

    %mode strength. It is the element that varies between open and rigid spheres
    %bn = Jn-dJn./dHn.*Hn;
    %alternative more simplified formulation using Wronskian relationship
    %(Teutsch's thesis)
%     bn = -2*1i./(dHn.*pi.*x);    
%     
%     p = p + (-1i).^n*exp(1i*n*(phi-phis))*bn;
    
    
    %Kaisser's expression (clarifications to be confirmed, not working)
%     if n<0
%         deltan = sin(phi*n);
%         deltans = sin(phis*n);
%     else
%         deltan = cos(phi*n);
%         deltans = cos(phis*n);
%     end
    %sine and cosine are swapped from Kaisser's definition as the
    %coordiante system is different (test)
%     if n<0
%         deltan = cos(phi*n);
%         deltans = cos(phis*n);
%     else
%         deltan = sin(phi*n);
%         deltans = sin(phis*n);
%     end
% 
%     
%     for imic=1:M
%         phin = sqrt((2-deltan(imic))/(2*pi));
%         phins = sqrt((2-deltans(imic))/(2*pi));
%         phin = deltan(imic);
%         phins = deltans(imic);
%         p(imic,:) = p(imic,:) + 2*pi*(1i)^(n+1).*phins.*phin.*exp(1i.*k.*sin(thetas).*z(imic))./ ...
%             (x.*cos(thetas).*dH2n);
%     end
    
%     phins = sqrt((2-deltans)/(2*pi));
%     phin = sqrt((2-deltan)/(2*pi));
%       %different coordinate system used so modifications are applied (e.g.
%         %cos is replaced by sin and vice versa
%     p = p + 2*pi*(1i)^(n+1).*phins.*phin.*exp(1i.*k.*sin(thetas).*z)./ ...
%     (x.*cos(thetas).*dH2n);

    % Teutsch's implementation based on eq (2.92) and asssuming infinitely
    % long cylinder, i.e. replacing integral with another summation
    
%     for m=-N:N %second integral for the cylinder's length
%         im=m+N+1;
%         km = k.*sin(theta(imic));   %I think it is wrong is not dependent on m
%         krho = k.*cos(theta(imic));
%         k_rhoxr = krho.*r_hat_norm2;
%         Hnkrhor = besselh(n,1,k_rhoxr);
%         dHnkrhor = 1/2.*(besselh(n-1,k_rhoxr) - besselh(n+1,k_rhoxr));
% 
%         pinc(in,1,imic,:) = (-1i)^n.*Jn.*exp(1i*n*(phi(imic)-phis(imic)));
%         pscat(in,im,imic,:) = (-1i)^n.*dJn.*exp(1i*n*(phi(imic)-phis(imic))) .* ...
%             Hnkrhor.*dHnkrhor.*exp(1i.*km.*z(imic));
% 
%     end


    % My attempt from Teustch derivation but for a plane wave impinging
    % from any angle and at any position on the rigid cylinder
 
    %the pressure does not use theta_P but it uses z, therefore inherently
    %cosidering the height at which the microphones are placed.
    
    dHnz = 1/2.*(besselh(n-1,2,krsinthetas) - besselh(n+1,2,krsinthetas));
    % effect of rygid cylinder bn using the Wronskian relation
    bn = 2./(1i.*pi.*krsinthetas.*dHnz);
    %ptot(:,:,in) = ((1i).^n.*exp(1i.*n.*(phi_P-phis_P))*bn).*e1;
    for imic=1:M       
        ptot(imic,:,in) = (1i).^n.*bn.*exp(1i.*n.*(phi_P(imic)-phis_P)).*e1(imic);
    end
    
end
% end
% pinc = squeeze(sum(pinc,1));
% pscat = squeeze(sum(sum(pscat,1),2));
% p = pinc + pscat;

p = squeeze(sum(ptot,3));

end