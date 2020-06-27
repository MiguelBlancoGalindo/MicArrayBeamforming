function [TR,W]=TriQuad(TR,W)
% Subdivide a triangular surface mesh using generalized triangular 
% quadrisection. Triangular quadrisection is a linear subdivision procedure
% which inserts new vertices into the input mesh at the edge midpoints 
% thereby producing four new faces for every face of the original mesh. 
% Illustration of this operation is provided below:
% 
%                     x3                        x3
%                    /  \      subdivision     /  \
%                   /    \         -->        v3__v2
%                  /      \                  / \  / \
%                x1________x2              x1___v1___x2
%
%                   Original vertices:    x1, x2, x3
%                   New vertices:         v1, v2, v3
% 
% In case of generalized triangular quadrisection, positions of the newly
% inserted vertices do not necessarily have to correspond to edge 
% midpoints, and could be varied by assigning weights to the vertices of 
% the original mesh. For example, let xi and xj be two vertices connected 
% by an edge, and suppose that Wi and Wj are the corresponding vertex 
% weights. The position of the new point on the edge (xi,xj) would be 
% defined as (Wi*xi+Wj*xj)/(Wi+Wj). Note that in order to avoid 
% degeneracies and self-intersections, all weights must be real numbers 
% greater than zero.
%
% INPUT ARGUMENTS:
%   - TR   : input mesh. TR can be specified as a TriRep object, structure
%            with fields 'faces' and 'vertices', or a cell so that 
%            TR={Tri X}, where X is a list of vertex co-ordinates and Tri 
%            is a list of faces.
%   - W    : optional input argument. N-by-1 array of NON-ZERO, POSITIVE 
%            vertex weights used during interpolation of the new vertices, 
%            where N is the total number of the original mesh vertices. 
%
% OUTPUT:
%   - TR  : subdivided mesh. Same format as the input.
%   - W   : new set of vertex weights.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: May.2012
%

% Get the list of vertex co-ordinates and list of faces
fmt=1;
if isa(TR,'TriRep')
    Tri=TR.Triangulation;
    X=TR.X;
elseif isstruct(TR) && sum(isfield(TR,{'vertices' 'faces'}))==2
    Tri=TR.faces;
    X=TR.vertices;
    fmt=2;
elseif iscell(TR) && numel(TR)==2
    Tri=TR{1};
    X=TR{2};
    fmt=3;
else
    error('Unrecognized input format')
end

% Make sure that the mesh is composed entirely of triangles 
if size(Tri,2)~=3 || ~isequal(round(Tri),Tri)
        error('Invalid entry for the list of faces')
end

% Check vertex weights
flag=false;
if nargin<2 || isempty(W), flag=true; end

if ~flag && (numel(W)~=size(X,1) || sum(W<=eps))
    error('W must be a N-by-1 array with non-negative entries, where N is the # of mesh vertices')
end
if ~flag, W=W(:)+eps; end

% Compute new vertex positions
if ~flag
    w=bsxfun(@rdivide,[W(Tri(:,1)),W(Tri(:,2))],W(Tri(:,1))+W(Tri(:,2)));
    V1=bsxfun(@times,X(Tri(:,1),:),w(:,1))+bsxfun(@times,X(Tri(:,2),:),w(:,2));
    w1=W(Tri(:,1)).*w(:,1)+W(Tri(:,2)).*w(:,2);
    
    w=bsxfun(@rdivide,[W(Tri(:,2)),W(Tri(:,3))],W(Tri(:,2))+W(Tri(:,3)));
    V2=bsxfun(@times,X(Tri(:,2),:),w(:,1))+bsxfun(@times,X(Tri(:,3),:),w(:,2));
    w2=W(Tri(:,2)).*w(:,1)+W(Tri(:,3)).*w(:,2);
    
    w=bsxfun(@rdivide,[W(Tri(:,3)),W(Tri(:,1))],W(Tri(:,3))+W(Tri(:,1)));
    V3=bsxfun(@times,X(Tri(:,3),:),w(:,1))+bsxfun(@times,X(Tri(:,1),:),w(:,2));
    w3=W(Tri(:,3)).*w(:,1)+W(Tri(:,1)).*w(:,2);
    
    W_new=[w1;w2;w3];
else
    V1=(X(Tri(:,1),:)+X(Tri(:,2),:))/2;
    V2=(X(Tri(:,2),:)+X(Tri(:,3),:))/2;
    V3=(X(Tri(:,3),:)+X(Tri(:,1),:))/2;
end
V=[V1;V2;V3];

% Remove repeating vertices 
[V,idx_unq,idx]=unique(V,'rows','stable'); % setOrder='stable' ensures that identical results (in terms of face connectivity) will be obtained for meshes with same topology
if ~flag
    W=[W;W_new(idx_unq)]; 
end

% Assign indices to the new triangle vertices
Nx=size(X,1);   % # of vertices
Nt=size(Tri,1); % # of faces

V1= Nx + idx(1:Nt);
V2= Nx + idx((Nt+1):2*Nt);
V3= Nx + idx((2*Nt+1):3*Nt);
clear idx

% Define new faces
T1= [Tri(:,1) V1 V3];
T2= [Tri(:,2) V2 V1];
T3= [Tri(:,3) V3 V2];
T4= [V1       V2 V3];
clear V1 V2 V3

T1=permute(T1,[3 1 2]);
T2=permute(T2,[3 1 2]);
T3=permute(T3,[3 1 2]);
T4=permute(T4,[3 1 2]);

Tri=cat(1,T1,T2,T3,T4);
Tri=reshape(Tri,[],3,1);

% New mesh
X=[X;V]; 
if fmt==1
    TR=TriRep(Tri,X);
elseif fmt==2
    TR.faces=Tri;
    TR.vertices=X;
else
    TR={Tri,X};
end
if nargout>1 && flag, W=ones(size(X,1),1); end

