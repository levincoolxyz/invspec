function varargout = lapbel(v,f)
% function [M,L] = lapbel(v,f)
% computes the Laplace-Beltrami operator of discrete surfaces
% 
% INPUTS
% v,f  - face-vertex data of the input surface
% 
% OUTPUTS
% M    - diagonal part (mass matrix)
% L    - symmetric negative definite part
%        thus M^(-1)*L is the Laplace-Beltrami operator
% 

numv = size(v,1); % number of vertices
numf = size(f,1); % number of faces
sparselimit = 5e2+1; % above which switch to sparse matrices
%% cot of angle at x and dual area associated to x in each face 
cotan = @(x,y,z) (y-x)'*(z-x)./norm(cross(y-x,z-x),2);
% dualA = @(x,y,z) (((y-x)'*(y-x))*cotan(z,x,y) + ((z-x)'*(z-x))*cotan(y,z,x))/8;
dualA = @(x,y,z) norm(cross(y-x,z-x),2)/6; % barycentric
%% construction (naive, way too many calls to @cotan)
if numv > sparselimit
  M = spalloc(numv,numv,numv);
else
  M = zeros(numv);
end
if (nargout > 1)
  if numv > sparselimit
    L = spalloc(numv,numv,7*numv);
  else
    L = zeros(numv);
  end
  for fi=1:numf
    for idx = 0:2
      i = f(fi,idx+1);
      j = f(fi,mod(idx+1,3)+1);
      k = f(fi,mod(idx+2,3)+1);
      L(i,k) = L(i,k) + .5*cotan(v(j,:)',v(k,:)',v(i,:)');
      L(i,i) = L(i,i) - .5*cotan(v(j,:)',v(k,:)',v(i,:)');
      L(i,j) = L(i,j) + .5*cotan(v(k,:)',v(i,:)',v(j,:)');
      L(i,i) = L(i,i) - .5*cotan(v(k,:)',v(i,:)',v(j,:)');
      M(i,i) = M(i,i) + dualA(v(i,:)',v(j,:)',v(k,:)');
    end
  end
  varargout{1} = M;
  varargout{2} = L;
else
  for fi=1:numf
    for idx = 0:2
      i = f(fi,idx+1);
      j = f(fi,mod(idx+1,3)+1);
      k = f(fi,mod(idx+2,3)+1);
      M(i,i) = M(i,i) + dualA(v(i,:)',v(j,:)',v(k,:)');
    end
  end
  varargout{1} = M;
end