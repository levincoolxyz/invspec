function [s,v] = meancurvflow(v0,f0,L0,M0,h)
% close all;
numv = size(v0,1);
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
update = @(v,vold) norm(vnorm(v - vold));

zentrum0 = volcenter(v0,f0);
v0 = v0 - repmat(zentrum0,numv,1);
% scl0 = calcvol(v0,f0)^(1/3);
scl0 = makeunitarea(v0,f0);

M = M0;
L = L0;
vold = v0;
v = ones(size(v0));
I = eye(size(v0,1));

% figure();trimesh(TriRep(f0,v0)); axis equal;

iter = 1;
while true
  for dim = 1:3
    v(:,dim) = (I - h*(M\L))\vold(:,dim);
  end
  v = v - repmat(volcenter(v,f0),numv,1);
%   scl = calcvol(v,f0)^(1/3);
  scl = makeunitarea(v,f0);
  
  v = v*scl/scl0;
%   figure();trimesh(TriRep(f0,v)); axis equal;
%   [M,L] = lapbel(v,f0); % traditional
  M = lapbel(v,f0); % conformal
%   [~,L] = lapbel(v,f0); % authalic
  dv = update(v,vold);
  fprintf('flow iter#%d; J = %g\n',iter,dv);
  iter = iter + 1;
  if dv <= 1e-1
    break;
  end
  vold = v;
end

s = diag(inv(M)*M0);
end
function [vol] = calcvol(v,f)
vol = 0;
for fi = 1:size(f,1)
  vi = f(fi,:);
  vol = vol + dot(v(vi(1),:),cross(v(vi(2),:),v(vi(3),:)))/6;
end
end
function scl = makeunitarea(v,f)
area = 0;
for fi = 1:size(f,1)
  vi = f(fi,:);
  area = area + norm(cross(v(vi(2),:) - v(vi(1),:),...
    v(vi(3),:) - v(vi(1),:)))/2;
end
scl = 1./sqrt(area);
end
function [zentrum,vol] = volcenter(v,f)
zentrum = zeros(1,3);
vol = 0;
for fi = 1:size(f,1)
  v1 = v(f(fi,1),:);
  v2 = v(f(fi,2),:);
  v3 = v(f(fi,3),:);
  vol = vol + det([v1;v2;v3]);
  zentrum = zentrum + (v1 + v2 + v3)*det([v1;v2;v3])/4;
end
zentrum = zentrum/vol;
end