function [OffsetPos, OffsetDistance] = ArrayOffset(OrigPos,sigma,ManifoldType)
% function [OffsetPos, OffsetDistance] = ArrayOffset(OrigPos,sigma,ManifoldType)
% function to calculate microphone array positioning offset 
% input parameters:
%   OrigPos: original microphone positions.
%   sigma: standard deviation of the offset. 
%   ManifoldType: type of array manifold. This will determine whether the 
%       offset set for each x, y, z coordinates (if 'freefield'), or 
%       contrained to the surface of a sphere ('spherical_scatterer') or 
%       cylinder ('cylindrical_scatterer').
%
% output arguments: 
%   OffsetPos: microphone array positions with offset
%   OffsetDistance: equivalent distance of offset introduced at each
%       microphone.

M = size(OrigPos,1);
ND = size(OrigPos,2);
Offset = zeros(M,ND);
if strcmp(ManifoldType, 'freefield')
    Offset = sigma*randn(size(OrigPos));
    OffsetDistance = sqrt(Offset(:,1).^2 + Offset(:,2).^2 + Offset(:,3).^2);
    OffsetPos = OrigPos + Offset;
elseif strcmp(ManifoldType, 'spherical_scatterer')
    OrigPosSph = zeros(M,ND);
    [OrigPosSph(:,1),OrigPosSph(:,2),OrigPosSph(:,3)] = car2sph(OrigPos(:,1),OrigPos(:,2),OrigPos(:,3));
    r = mean(OrigPosSph(:,3));
    %conversion from Cartesian distance std to angle offset in Sph coord
    %including conversion from 3 components(x,y,z) to 2(phi,theta) in offset.
    sigmaSph = sqrt(3/2).*sigma/r.*180./pi;
    OffsetSph = zeros(M,ND);
    OffsetSph(:,1) = sigmaSph*randn(M,1);
    OffsetSph(:,2) = sigmaSph*randn(M,1);
    %compensation for smaller apparent radius near the poles, i.e. azimuth
    %offset should be larger near the poles for a fixed distance offset. 
    OffsetSph(:,1) = OffsetSph(:,1)./cosd(OrigPosSph(:,2));
    OffsetPosSph = OrigPosSph + OffsetSph;
    [OffsetPos(:,1),OffsetPos(:,2),OffsetPos(:,3)] = sph2car(OffsetPosSph(:,1),OffsetPosSph(:,2),OffsetPosSph(:,3));
    [OffsetDistance,OffsetDistanceLong,OffsetDistanceLat] = spherror(OrigPosSph(:,1),OrigPosSph(:,2),OffsetPosSph(:,1),OffsetPosSph(:,2));
    OffsetDistance = OffsetDistance.*pi/180.*r;
    
elseif strcmp(ManifoldType, 'cylindrical_scatterer')
    OrigPosCyl = zeros(M,ND);
    [OrigPosCyl(:,1),OrigPosCyl(:,2),OrigPosCyl(:,3)] = car2cyl(OrigPos(:,1),OrigPos(:,2),OrigPos(:,3));
    rho = mean(OrigPosCyl(:,3));
    %conversion from Cartesian distance std to angle offset in Cyl coord
    %including conversion from 3 components(x,y,z) to 2(phi,theta) in offset.
    sigmaCyl = sqrt(3/2).*sigma;
    OffsetCyl = zeros(M,ND);
    OffsetCyl(:,1) = repmat(sigmaCyl.*1/rho.*180./pi,M,1).*randn(M,1);
    OffsetCyl(:,2) = repmat(sigmaCyl,M,1).*randn(M,1);
    OffsetPosCyl = OrigPosCyl + OffsetCyl;
    [OffsetPos(:,1),OffsetPos(:,2),OffsetPos(:,3)] = cyl2car(OffsetPosCyl(:,1),OffsetPosCyl(:,2),OffsetPosCyl(:,3));
    OffsetCylAzDis = OffsetCyl(:,1).*pi/180.*rho;
    OffsetDistance = sqrt(OffsetCylAzDis.^2 + OffsetCyl(:,2).^2 + OffsetCyl(:,3).^2);
end

