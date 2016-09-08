function [s,v,M] = meancurvflow(v0,f0,h,type,L0,M0,imax,gif)
% function [s,v,M] = meancurvflow(v0,f0,h,type,L0,M0,imax,gif)
% Mean Curvature Flow of discrete surfaces
% 
% INPUTS
% v0,f0    - face-vertex data of the initial surface
% h        - step size of the flow
% type     - type of flow wanted ('t','c','a')
% L0,M0    - (optional) precomupted initial Laplace-Beltrami operator
% imax     - (optional) iteration limit
% gif      - flag 1 to produce animated results, default to 0 otherwise
% 
% OUTPUTS
% s        - resultant conformal factors
% v        - resultant vertex coordinates
% 

if (nargin<3) || isempty(h), h = 1; end
if (nargin<4) || isempty(type), type = 'c'; end % default to cMCF
if (nargin<5) || isempty(L0), [M0,L0] = lapbel(v0,f0); end
if (nargin<6) || isempty(M0), M0 = lapbel(v0,f0); end
if (nargin<7) || isempty(imax), imax = 1e3; end
if (nargin<8) || isempty(gif), gif = 0; end

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
% I = eye(size(v0,1));

if gif
  figure();trimesh(f0,v0(:,1),v0(:,2),v0(:,3)); axis equal;
  saveas(gcf,'mcf000.png'); close(gcf);
end

iter = 1;
sphericity_old = std(vnorm(v0));
while iter<=imax
  if numv > 500
    M = sparse(M);
    L = sparse(L);
  end
  for dim = 1:3
%     A = (I - h*(M\L));
%     b = vold(:,dim);
    A = (M - h*L);
    b = (M*vold(:,dim));
    v(:,dim) = A\b;
%     l = chol(A,'lower');
%     w = l\b;
%     v(:,dim) = l'\w;
  end
  v = v - repmat(volCenter(v,f0),numv,1);
%   scl = calcVol(v,f0)^(1/3);
  scl = makeUnitArea(v,f0);
  v = v*scl/scl0;
  

if gif
  figure();trimesh(f0,v(:,1),v(:,2),v(:,3)); axis equal;
  saveas(gcf,num2str(iter,'mcf%03d.png')); close(gcf);
end

  if type == 't'
    [M,L] = lapbel(v,f0); % traditional MCF
  elseif type == 'c'
    M = lapbel(v,f0); % conformal MCF
  elseif type == 'a'
    [~,L] = lapbel(v,f0); % authalic MCF
  else
    error('wat?');
  end
  
  sphericity = std(vnorm(v));
  dv = vdiff(v,vold);
  fprintf('MCF iter#%d: |dv| = %g; sphericity = %g\n',iter,dv, sphericity);
  iter = iter + 1;
  if abs(sphericity_old - sphericity) <= 2e-5 && sphericity <= 2e-2
    break;
  end
  if sphericity > 2
    warning('unstable mean curvature flow');
    break;
  end
  vold = v;
  sphericity_old = sphericity;
end

% s = 1./diag(M\M0);
% s = diag(M0\M);
s = diag(M)./diag(M0);

s = full(s);

if gif
  unix(['convert -delay 10 -loop 0 mcf*.png mcf/mcf' num2str(h,'%g') '.gif']);
  unix('rm mcf*.png');
end

% figure();trimesh(f0,v(:,1),v(:,2),v(:,3)); axis equal;
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