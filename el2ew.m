function [M,L,trineq] = el2ew(f,elsq)
% function [M,L] = el2ew(f,el)
% Get edge weights (cotan weights, massless laplacian) from edge lengths (metric)
% 
% INPUT
% f      - (numf x 3) face connectivity
% elsq   - (nume x 1) edge lengths squared in isedge order (lowertriangular)
% OUTPUT
% M      - (numv x numv) (diagonal) mass matrix
% L      - (numv x numv) cotan weight matrix
% trineq - (numf x 1) check which face defies triangle inequality
% 

numv = numel(unique(f(:)));
numf = size(f,1); % number of faces
sparselimit = 1e3; % above which switch to sparse matrices
%% fuss with (is)edge indices
isedge = zeros(numv);
for fi = 1:numf
  for idx = 0:2
    i = f(fi,idx+1);
    j = f(fi,mod(idx+1,3)+1);
    isedge(i,j) = 1;
  end
end
isedge = tril(isedge); % reduce redundancy
isedge = find(isedge); % linear indices
eidx1 = @(i,j) sub2ind([numv numv],i,j)*(i>j) + sub2ind([numv numv],j,i)*(j>i);
eidx = @(i,j) find(isedge == eidx1(i,j));
%% compute cotan weights and (barycentric) masses
% cotan at vertex between edge x and edge y
cotan = @(x,y,z) (x + y - z)/sqrt(-x^2 + 2*x*y + 2*x*z - y^2 + 2*y*z - z^2);
pma = @(x,y,z,a) sum([x y z])/2 - a; % perimeter minus a
el = sqrt(elsq); % get actual edge length for area
dualA = @(x,y,z) sqrt(pma(x,y,z,0)*pma(x,y,z,x)*pma(x,y,z,y)*pma(x,y,z,z))/3; % barycentric
if numv > sparselimit
  L = spalloc(numv,numv,7*numv);
  M = spalloc(numv,numv,numv);
else
  L = zeros(numv);
  M = zeros(numv);
end
if nargout > 2
  trineq = zeros(size(f,1),1);
end
for fi=1:numf
  for idx = 0:2
    i = f(fi,idx+1);
    j = f(fi,mod(idx+1,3)+1);
    k = f(fi,mod(idx+2,3)+1);
    I = eidx(j,k);
    J = eidx(i,k);
    K = eidx(i,j);
    L(i,k) = L(i,k) + .5*cotan(elsq(I),elsq(K),elsq(J));
    L(i,i) = L(i,i) - .5*cotan(elsq(I),elsq(K),elsq(J));
    L(i,j) = L(i,j) + .5*cotan(elsq(J),elsq(I),elsq(K));
    L(i,i) = L(i,i) - .5*cotan(elsq(J),elsq(I),elsq(K));
    M(i,i) = M(i,i) + dualA(el(J),el(I),el(K));
    
    if nargout > 2
      % check triangle inequality
      elf = sort([el(I),el(J),el(K)]);
      trineq(fi) = elf(1) + elf(2) <= elf(3);
    end
  end
end
end