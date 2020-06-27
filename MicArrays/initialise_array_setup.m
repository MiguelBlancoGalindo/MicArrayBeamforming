function Config = initialise_array_setup(Config)
% function Config = initialise_array_setup(Config)
% function that calculates the sensor spacing, radius and/or number of
% sensors required to satisfy one of the three criteria in
% "Config.Array.Criteria" with the circular array as the reference

% input and output parameter:
%   Config: configuration structure containing all the settings.

Config.Array.Spacing=[];
Config.Array.Orientation='xy';
if strcmp(Config.Array.Geo,'stacked_circular') || strcmp(Config.Array.Geo,'stacked_circular_rigid_cylinder')
    if ~isfield(Config.Array,'NCircles')
        Config.Array.NCircles = 2;
    end
    if ~isfield(Config.Array,'StackedCircleSpacing')
        Config.Array.StackedCircleSpacing = 0.05;
    end
    if ~isfield(Config.Array,'MultiCircleAlignment')
        Config.Array.MultiCircleAlignment = 'aligned';
    end
elseif strcmp(Config.Array.Geo,'multicircular')
    Config.Array.NCircles = 4;
    Config.Array.MultiCircleAlignment = 'notaligned';
else
    Config.Array.NCircles = [];
    Config.Array.MultiCircleAlignment = [];
end

if strcmp(Config.Array.Criteria, 'spacing')
    %the spacing of the geometries is kept fixed
    if strcmp(Config.Array.Geo,'linear') || strcmp(Config.Array.Geo,'linear_bs')
        Config.Array.Spacing=2*Config.Array.r*sin(pi/Config.Array.M);
        Config.Array.r=[];
        Config.Array.Orientation = 'x';
    elseif strcmp(Config.Array.Geo,'linear_ef')
        Config.Array.Spacing=2*Config.Array.r*sin(pi/Config.Array.M);
        Config.Array.r=[];
        Config.Array.Orientation = 'y';
    elseif strcmp(Config.Array.Geo,'rectangular') || strcmp(Config.Array.Geo,'square')
        Config.Array.Spacing=2*Config.Array.r*sin(pi/Config.Array.M);
        Config.Array.r=[];
    elseif strcmp(Config.Array.Geo,'dualcircular')
        r1=Config.Array.r*sin(pi/Config.Array.M)/sin(pi/(Config.Array.M/2));
        Config.Array.r = [r1 r1/0.8];    %may be changed
    elseif strcmp(Config.Array.Geo,'spherical') || strcmp(Config.Array.Geo,'spherical_rigid_sphere')
        if Config.Array.M==48
            Config.Array.r=0.055;    %calculated by trying different Config.Array.r to give the same spacing as that of the dual-circular
        elseif Config.Array.M==12
        end
    end
elseif strcmp(Config.Array.Criteria, 'aperture')
    %the aperture of the geometries is kept fixed
    if strcmp(Config.Array.Geo,'linear') || strcmp(Config.Array.Geo,'linear_bs')
        Config.Array.Spacing=2*Config.Array.r/(Config.Array.M-1);
        Config.Array.r=[];
        Config.Array.Orientation = 'x';
    elseif strcmp(Config.Array.Geo,'linear_ef')
        Config.Array.Spacing=2*Config.Array.r/(Config.Array.M-1);
        Config.Array.r=[];
        Config.Array.Orientation = 'y';
    elseif strcmp(Config.Array.Geo,'rectangular')
        if sqrt(Config.Array.M)~=floor(sqrt(Config.Array.M))
            My=1;
            Mx=Config.Array.M;
            exact = zeros(Config.Array.M,1);
            while My<Mx
                Mx=floor(Config.Array.M/My);
                if mod(Config.Array.M,My)==0
                    exact(My)=My;
                end
                My=My+1;
            end
            My=max(exact);
            Mx=Config.Array.M/My;
        else
            My=sqrt(Config.Array.M);
            Mx=My;
        end
        Config.Array.Spacing=2*Config.Array.r/(Mx-1);
        Config.Array.r=[];
    elseif strcmp(Config.Array.Geo,'square')
        Mx = sqrt(Config.Array.M);
        if round(Mx.^2)==Config.Array.M
            Config.Array.Spacing = 2*Config.Array.r/(Mx-1);
        else
            error('the number of sensor is not an even square');
        end
    elseif strcmp(Config.Array.Geo,'dualcircular')
        Config.Array.r=[0.8*Config.Array.r Config.Array.r];
    end
elseif strcmp(Config.Array.Criteria, 'aperture_and_spacing')
    %the aperture and the spacing are kept fixed by varying Config.Array.M
    if strcmp(Config.Array.Geo,'linear') || strcmp(Config.Array.Geo,'linear_bs')
        Config.Array.Spacing=2*Config.Array.r*sin(pi/(Config.Array.M));
        Config.Array.M=floor(2*Config.Array.r/Config.Array.Spacing) + 1;
        Config.Array.Orientation = 'x';
    elseif strcmp(Config.Array.Geo,'linear_ef')
        Config.Array.Spacing=2*Config.Array.r*sin(pi/(Config.Array.M));
        Config.Array.M=floor(2*Config.Array.r/Config.Array.Spacing) + 1;
        Config.Array.Orientation = 'y';
    elseif strcmp(Config.Array.Geo,'planar')
        Config.Array.Spacing=2*Config.Array.r*sin(pi/(Config.Array.M));
        l=sqrt(2)*Config.Array.r;
        M_side = floor(l/Config.Array.Spacing) + 1;
        Config.Array.M=M_side^2;
        Config.Array.r=[];
    elseif strcmp(Config.Array.Geo,'dualcircular')
        Config.Array.M = 2*Config.Array.M;
        Config.Array.r=[0.8*Config.Array.r Config.Array.r];
    elseif strcmp(Config.Array.Geo,'spherical') || strcmp(Config.Array.Geo,'spherical_rigid_sphere')
        Config.Array.M=4*Config.Array.M;  %rule of thumb so far
    end
end    