function [fv,W]=QuadQuad(fv,W)
% Subdivide a quadrilateral surface mesh using generalized quadrilateral 
% quadrisection. This procedure is similar to generalized triangular 
% quadrisection described in the documentation of 'TriQuad' function, the 
% only difference being that it uses quad surface meshes as input.
%
% INPUT ARGUMENTS:
%   - fv  :  input mesh specified as a face-vertex structure with fields 
%            'vertices' and 'faces' so that fv.vertices contains a N-by-3 
%            list of vertex co-ordinates and fv.faces contains a M-by-4 
%            list of faces. Alternatively, fv can be specified as cell so 
%            that fv={F X}, where X is a N-by-3 list of vertex co-ordinates 
%            and F is a M-by-4 list of faces.
%   - W    : optional input argument. N-by-1 array of NON-ZERO, POSITIVE 
%            vertex weights used during interpolation of the new vertices,
%            where N is the total number of the original mesh vertices. 
%
% OUTPUT:
%   - fv  : subdivided mesh.
%   - W   : new set of vertex weights.   
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: March, 2015
%

% Get the list of vertex co-ordinates and list of faces
fmt=true;
if isstruct(fv) && sum(isfield(fv,{'vertices' 'faces'}))==2
    F=fv.faces;
    X=fv.vertices;
elseif iscell(fv) && numel(fv)==2
    F=fv{1};
    X=fv{2};      
    fmt=false;
else
    error('Unrecognized input format')
end

% Make sure that the mesh is composed entirely of quads 
if size(F,2)~=4 || ~isequal(round(F),F)
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
    w=bsxfun(@rdivide,[W(F(:,1)),W(F(:,2))],W(F(:,1))+W(F(:,2)));
    V1=bsxfun(@times,X(F(:,1),:),w(:,1))+bsxfun(@times,X(F(:,2),:),w(:,2));
    w1=W(F(:,1)).*w(:,1)+W(F(:,2)).*w(:,2);
    
    w=bsxfun(@rdivide,[W(F(:,2)),W(F(:,3))],W(F(:,2))+W(F(:,3)));
    V2=bsxfun(@times,X(F(:,2),:),w(:,1))+bsxfun(@times,X(F(:,3),:),w(:,2));
    w2=W(F(:,2)).*w(:,1)+W(F(:,3)).*w(:,2);
    
    w=bsxfun(@rdivide,[W(F(:,3)),W(F(:,4))],W(F(:,3))+W(F(:,4)));
    V3=bsxfun(@times,X(F(:,3),:),w(:,1))+bsxfun(@times,X(F(:,4),:),w(:,2));
    w3=W(F(:,3)).*w(:,1)+W(F(:,4)).*w(:,2);
    
    w=bsxfun(@rdivide,[W(F(:,4)),W(F(:,1))],W(F(:,4))+W(F(:,1)));
    V4=bsxfun(@times,X(F(:,4),:),w(:,1))+bsxfun(@times,X(F(:,1),:),w(:,2));
    w4=W(F(:,4)).*w(:,1)+W(F(:,1)).*w(:,2);
    
    w=bsxfun(@rdivide,[W(F(:,1)) W(F(:,2)) W(F(:,3)) W(F(:,4))],...
                      W(F(:,1))+W(F(:,2)) + W(F(:,3))+W(F(:,4)));
    V5=bsxfun(@times,X(F(:,1),:),w(:,1))+bsxfun(@times,X(F(:,2),:),w(:,2))+...
       bsxfun(@times,X(F(:,3),:),w(:,3))+bsxfun(@times,X(F(:,4),:),w(:,4));
    w5=W(F(:,1)).*w(:,1)+W(F(:,2)).*w(:,2)+W(F(:,3)).*w(:,3)+W(F(:,4)).*w(:,4);
   
    W_new=[w1;w2;w3;w4;w5];
else
    V1=(X(F(:,1),:)+X(F(:,2),:))/2;
    V2=(X(F(:,2),:)+X(F(:,3),:))/2;
    V3=(X(F(:,3),:)+X(F(:,4),:))/2;
    V4=(X(F(:,4),:)+X(F(:,1),:))/2;
    V5=(V1+V3)/2;
end
V=[V1;V2;V3;V4;V5];

% Remove repeating vertices 
[V,idx_unq,idx]=unique(V,'rows','stable'); % setOrder='stable' ensures that identical results (in terms of face connectivity) will be obtained for meshes with same topology
if ~flag
    W=[W;W_new(idx_unq)]; 
end

% Assign indices to the new vertices
Nx=size(X,1);   % # of vertices
Nt=size(F,1); % # of faces

V1= Nx + idx(1:Nt);
V2= Nx + idx((Nt+1):2*Nt);
V3= Nx + idx((2*Nt+1):3*Nt);
V4= Nx + idx((3*Nt+1):4*Nt);
V5= Nx + idx((4*Nt+1):5*Nt);
clear idx

% Define new faces
T1= [F(:,1) V1 V5 V4];
T2= [V1 F(:,2) V2 V5];
T3= [V5 V2 F(:,3) V3];
T4= [V4 V5 V3 F(:,4)];
clear V1 V2 V3

T1=permute(T1,[3 1 2]);
T2=permute(T2,[3 1 2]);
T3=permute(T3,[3 1 2]);
T4=permute(T4,[3 1 2]);

F=cat(1,T1,T2,T3,T4);
F=reshape(F,[],4,1);

% New mesh
clear fv
if fmt
    fv.faces=F;
    fv.vertices=[X;V];
else
    fv={F [X;V]};
end
if nargout>1 && flag, W=ones(size(X,1),1); end

