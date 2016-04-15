function [s,v] = meancurvflow(v0,f0,L0,M0,h)
close all;

vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
update = @(v,vold) norm(vnorm(v - vold));
scl = @(M) (sum(diag(M)));
scl0 = scl(M0);
M = M0;
L = L0;
vold = v0;
v = ones(size(v0));
I = eye(size(v0,1));
iter = 1;
while true
  for dim = 1:3
    v(:,dim) = (I - h*(M\L))\vold(:,dim);
  end
  v = v/scl(M)*scl0;
  figure();trimesh(TriRep(f0,v)); axis equal;
%   [M,L] = lapbel(v,f0); % traditional
  M = lapbel(v,f0); % conformal
%   [~,L] = lapbel(v,f0); % authalic
  dv = update(v,vold);
  fprintf('flow iter#%d; J = %g\n',iter,dv);
  iter = iter + 1;
  if dv <= 1e-5
    break;
  end
  vold = v;
end

s = diag(inv(M)*M0);