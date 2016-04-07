function [M,L] = lapbel(v,f)
numv = size(v,1); % number of vertices
numf = size(f,1); % number of faces
%% cot of angle at x and dual area associated to x in each face 
cotan = @(x,y,z) (y-x)'*(z-x)./norm(cross(y-x,z-x),2);
dualA = @(x,y,z) ((y-x)'*(y-x))*cotan(z,x,y) + ((z-x)'*(z-x))*cotan(y,z,x);
%% construction (naive, way too many calls to @cotan)
% L = spalloc(numv,numv,6*numv);
% M = spalloc(numv,numv,numv);
L = zeros(numv);
M = zeros(numv);
for fi=1:numf
  for idx = 0:2
    i = f(fi,idx+1);
    j = f(fi,mod(idx+1,3)+1);
    k = f(fi,mod(idx+2,3)+1);
    L(i,k) = L(i,k) + .5*cotan(v(j,:)',v(k,:)',v(i,:)');
    L(i,i) = L(i,i) - .5*cotan(v(j,:)',v(k,:)',v(i,:)');
    L(i,j) = L(i,j) + .5*cotan(v(k,:)',v(i,:)',v(j,:)');
    L(i,i) = L(i,i) - .5*cotan(v(k,:)',v(i,:)',v(j,:)');
    M(i,i) = M(i,i) + .125*dualA(v(i,:)',v(j,:)',v(k,:)');
  end
end
% L = sparse(L);
% M = sparse(M);