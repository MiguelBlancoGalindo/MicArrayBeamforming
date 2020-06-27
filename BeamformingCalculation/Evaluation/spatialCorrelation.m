%function that calculates the spatial correlation between theoretical and
%synthesised responses as per:
    %Moreau, Sébastien, Jérôme Daniel, and Stéphanie Bertet. 
    %"3D sound field recording with higher order ambisonics–Objective measurements 
    %and validation of a 4th order spherical microphone." 120th Convention of the AES. 2006

% d_d: [L x F x S] desired (or target) directivity response. L is the number of
% directions in space, F is the number of frequencies and S is the number
% of simulations
% d: [L x F x S] synthesised directivity response. Can also be a column
% vector
% additionally weights of curvature over the sphere can be passed with 'w' 
% flag followed by the weight vector

function R = spatialCorrelation(d_d, d,varargin)

nD = ndims(d);
L = size(d,1);
F = size(d,2);
S = size(d,3);
if numel(d_d)~=numel(d)
    if size(d_d,2)==1
        d_d = repmat(d_d,1,F,S);
    else
        if size(d_d,3)==1
            d_d = repmat(d_d,1,1,S);
        end
    end
   
end
dSize = size(shiftdim(d,1));
dSize(end) = 1;
R = zeros(dSize(1:nD));
nvarargin = length(varargin);
if nvarargin>0
    if strcmp(varargin{1},'w')
        w = varargin{2};
    end
else
    w = 1/L.*ones(L,1);
end

for ifreq=1:F
    for s=1:S
        num = innerProduct(d_d(:,ifreq,s),d(:,ifreq,s),w);
        R(ifreq,s) = num./(sqrt(innerProduct(d_d(:,ifreq,s),d_d(:,ifreq,s),w)).*sqrt(innerProduct(d(:,ifreq,s),d(:,ifreq,s),w)));
    end
end

function iProduct = innerProduct(x,y,w)
    iProduct = x'*diag(w)*y;
end


end
