function TR=SubdivideSphericalMesh(TR,k,opt)
% Subdivide triangular/quadrilateral mesh representing the surface of a 
% unit sphere k times using triangular/quadrilateral quadrisection (see
% function 'TriQuad'/'QuadQuad' for more info).
%
% INPUT ARGUMENTS:
%   - TR   : input mesh. TR must be specified as a TriRep object, a
%            structure with fields 'faces' and 'vertices', or a cell {F,V},
%            where F is a list faces and V is a list of vertices.
%   - k    : desired number of subdivisions. k=1 is default.
%   - opt  : optional input argument specifying how the newly inserted
%            vertices should be constrained to the spherical surface.
%            When opt=true, the newly inserted vertices are re-projected 
%            onto the surface of a unit sphere after every iteration. When
%            opt=false, projection onto the sphere is performed after the 
%            sequence of all k subdivisions had been completed. opt=true is
%            the default setting.
%
% OUTPUT:
%   - TR  : subdivided mesh. Same format as input.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: June.2012
%

if nargin<2 || isempty(k), k=1; end
if nargin<3 || isempty(opt), opt=true; end

% Get the data structure
fmt=1;
if isa(TR,'TriRep')
    F=TR.Triangulation;
    X=TR.X;
elseif isstruct(TR) && sum(isfield(TR,{'vertices' 'faces'}))==2
    F=TR.faces;
    X=TR.vertices;
    fmt=2;
elseif iscell(TR) && numel(TR)==2
    F=TR{1};
    X=TR{2};
    fmt=3;
else
    error('Unrecognized input format')
end

% Make sure that the input mesh is either triangular or quadrilateral 
if ~(size(F,2)==3 || size(F,2)==4) || ~isequal(round(F),F)
    error('''SubdivideSphericalMesh'' only works with triangular and quadrilateral meshes')
end

% If k<1 normalize the norm of vertex co-ordinates to unity and return
k=round(k(1));
if k<1 
    X_L2=sqrt(sum(X.^2,2));
    X=bsxfun(@rdivide,X,X_L2);
    if fmt==1
        TR=TriRep(F,X);
    elseif fmt==2
        TR.vertices=X;
    else
        TR{2}=X;
    end
    return
end


% Spherical subdivision
for i=1:k
    
    % Subdivide the mesh
    if size(F,2)==3
        TR=TriQuad({F X});
    else
        TR=QuadQuad({F X});
    end
    F=TR{1};
    X=TR{2};
        
    % Project the points onto the surface of the unit sphere
    if opt || i==k
        X_L2=sqrt(sum(X.^2,2));
        X=bsxfun(@rdivide,X,X_L2);
    end

end

clear TR
if fmt==1
    TR=TriRep(F,X);
elseif fmt==2
    TR.faces=F;
    TR.vertices=X;
else
    TR{1}=F;
    TR{2}=X;
end
