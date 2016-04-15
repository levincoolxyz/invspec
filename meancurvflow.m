function [s,M] = meancurvflow(v0,f0,L0,M0,h)

vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
cost = @(vnew,v) norm(vnorm(vnew - v));
M = M0;
L = L0;
vold = v0;
v = ones(size(v0));
I = eye(size(v0,1));
iter = 1;
while true
  for dim = 1:3
    v(:,dim) = (I - h*inv(M)*L)\vold(:,dim);
  end
%   trimesh(TriRep(f0,v));
%   [M,L] = lapbel(v,f0);
  M = lapbel(v,f0);
  J = cost(v,vold);
  fprintf('flow iter#%d; J = %g\n',iter,J);
  iter = iter + 1;
  if J <= 1e-5
    break;
  end
  vold = v;
end

s = diag(inv(M)*M0);