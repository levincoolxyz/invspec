function [s,v] = meancurvflow(v0,f0,h,type,L0,M0)
% function [s,v] = meancurvflow(v0,f0,h,type,L0,M0)
% Mean Curvature Flow of discrete surfaces
% 
% INPUTS
% v0,f0    - face-vertex data of the initial surface
% h        - step size of the flow
% L0,M0    - precomupted initial Laplace-Beltrami operator
% type     - type of flow wanted ('t','c','a')
% 
% OUTPUTS
% s        - resultant conformal factors
% v        - resultant vertex coordinates
% 

if (nargin<3) || isempty(h), h = 100; end
if (nargin<4) || isempty(type), type = 'c'; end % default to cMCF
if (nargin<5) || isempty(L0), [M0,L0] = lapbel(v0,f0); end

% close all;
numv = size(v0,1);
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
vdiff = @(v,vold) norm(vnorm(v - vold));

zentrum0 = volCenter(v0,f0);
v0 = v0 - repmat(zentrum0,numv,1);
% scl0 = calcVol(v0,f0)^(1/3);
scl0 = makeUnitArea(v0,f0);

M = M0;
L = L0;
vold = v0;
v = ones(size(v0));
I = eye(size(v0,1));

% figure();trimesh(TriRep(f0,v0)); axis equal;
sphericity_old = std(vnorm(v0));

iter = 1;
while true
  for dim = 1:3
    A = (I - h*(M\L));
    v(:,dim) = A\vold(:,dim);
  end
  v = v - repmat(volCenter(v,f0),numv,1);
%   scl = calcVol(v,f0)^(1/3);
  scl = makeUnitArea(v,f0);
  
  v = v*scl/scl0;
%   figure(); trimesh(TriRep(f0,v)); axis equal; pause(.5)
  sphericity = std(vnorm(v));

  if type == 't'
    [M,L] = lapbel(v,f0); % traditional
  elseif type == 'c'
    M = lapbel(v,f0); % conformal
  elseif type == 'a'
    [~,L] = lapbel(v,f0); % authalic
  else
    error('wat?');
  end
  
  dv = vdiff(v,vold);
  fprintf('flow iter#%d; v difference = %g\n',iter,dv);
  iter = iter + 1;
  if abs(sphericity_old - sphericity) <= 1e-5
    break;
  end
  vold = v;
  sphericity_old = sphericity;
end

% s = diag(inv(M)*M0);
s = diag(M\M0);
end

% function [vol] = calcVol(v,f)
% vol = 0;
% for fi = 1:size(f,1)
%   vi = f(fi,:);
%   vol = vol + dot(v(vi(1),:),cross(v(vi(2),:),v(vi(3),:)))/6;
% end
% end

function scl = makeUnitArea(v,f)
area = 0;
for fi = 1:size(f,1)
  vi = f(fi,:);
  area = area + norm(cross(v(vi(2),:) - v(vi(1),:),...
    v(vi(3),:) - v(vi(1),:)))/2;
end
scl = 1./sqrt(area);
end

function [zentrum,vol] = volCenter(v,f)
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