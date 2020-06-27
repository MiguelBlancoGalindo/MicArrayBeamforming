function fv=QuadCubeMesh
% Get a quad mesh of a cube whose vertices lie on the surface of a unit
% sphere.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: March, 2015  
%


X=[1 1; -1 1; -1 -1; 1 -1];
[Xd,Xu]=deal(X);
Xd(:,3)=-1;
Xu(:,3)=1;
X=[Xd;Xu];

F=[1 5 8 4; ...
   2 6 5 1; ...
   3 7 6 2; ...
   4 8 7 3; ...
   5 6 7 8; ...
   1 4 3 2];

fv.faces=F;
fv.vertices=X/sqrt(3);