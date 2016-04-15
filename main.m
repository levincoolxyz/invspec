clear;
%% initialization
input_case = 4; % 1 - octa; 2 - icosa; 3 - *.obj; 4 - sphere of ~ssize
target_case = 4; % 1 - octa; 2 - conf defms; 3 - *.obj; 4 - spharm defms
imax = 5e3; % gradient descent maximum iterations
aC = .5; bC = .8; etolC = 1e-3; % Conformal gradient descent control
aS = .5; bS = .4; etolS = 1e-3; % invSpec gradient descent control
numeig = 0; % number of eigenvalues used, 0 means full input
rng(1432543); % rand seed
purt = .6; % scaling coefficient used to control target purtabation
ssize = 200;
%% some spherical harmonics
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
Y33 = @(v) ((v(:,1).^2-3*v(:,2).^2).*v(:,1))./vnorm(v);
Y20 = @(v) (2*v(:,3).^2-v(:,1).^2-v(:,2).^2)./vnorm(v);
Y10 = @(v) v(:,3)./vnorm(v);
sphar = @(v) abs(Y33(v))*purt^3;
%% input mesh
% regular octahedron (1)
if input_case == 1
  v = [0 0 1; 0 1 0; 1 0 0 ;-1 0 0; 0 -1 0; 0 0 -1];
  f=fliplr(convhulln(v));
  
% regular icosahedron (2)
% cf. Anton Semechko (a.semechko@gmail.com)
elseif input_case == 2
  t=(1+sqrt(5))/2; % golden ratio
  v=[0 1 t];
  s=[1 1 1; 1 1 -1; 1 -1 -1; 1 -1 1];
  v=repmat(v,[4 1]).*s;
  v=[v;circshift(v,[0 -1]);circshift(v,[0 -2])];
  v_L2=sqrt(sum(v.^2,2));
  v=bsxfun(@rdivide,v,v_L2);
  clear s; clear t;
  f=fliplr(convhulln(v));
  
% import wavefront object file (3)
elseif input_case == 3
  filename = 'sphere_small';
  fid = fopen(['../meshes/' filename '.obj'],'rt');
  [v,f] = readwfobj(fid);

% sphere of chosen size (4)
elseif input_case == 4
%   [x,y,z]=sphere(floor(sqrt(ssize)));
%   v = unique([x(:) y(:) z(:)],'rows');
  v = ParticleSampleSphere('Vo',RandSampleSphere(ssize));
  f = fliplr(convhulln(v));
end

numv = size(v,1); % number of vertices
numf = size(f,1); % number of faces
numeig = numv*(numeig <= 0) + min(numv,numeig)*(numeig > 0);
%% when is there an edge (mild redundancy)
isedge = zeros(numv);
for fi = 1:numf
  for idx = 0:2
    i = f(fi,idx+1);
    j = f(fi,mod(idx+1,3)+1);
    isedge(i,j) = 1;
  end
end
isedge = triu(isedge); % reduce redundancy
isedge = find(isedge); % linear indices
%% compute initial edge lengths squared
elsq0 = zeros(numv);
for i = 1:numv
  for j = (i+1):numv % skipping the symmetric lower triangular part
    elsq0(i,j) = sum((v(i,:)-v(j,:)).^2);
  end
end
% elsq0 = elsq0.*isedge;
elsq0 = elsq0(isedge); % linear indices
%% compute laplacian
[M,L] = lapbel(v,f);
D_0 = eigvf(L,M,numeig);
%% compute mean curvature vertex normal
Hn = .5*[inv(M)*L*v(:,1) inv(M)*L*v(:,2) inv(M)*L*v(:,3)];
H = vnorm(Hn);
vn = Hn./repmat(H,1,3);
%% target spectrum
% skewed octahedron
if target_case == 1
  v_T = [0 0 .5; 0 1 0; 2 0 0 ;-2 0 0; 0 -3 0; 0 0 -.5];
  f_T = fliplr(convhulln(v));
  [M_T,L_T] = lapbel(v_T,f_T);
  D_T = eigvf(L_T,M_T,numeig);

% random small prescribed conformal deformation
elseif target_case == 2
  s_T = exp(-rand(numv,1)*purt);
  D_Tp = eigvf(L,diag(1./s_T)*M,numeig);
  f_T = f;
  conf_T = sqrt(kron(1./s_T',1./s_T));
%   elsq_T = elsq0.*conf_T;
  elsq_T = elsq0.*conf_T(isedge); % linear indices
  [Jc_Thist,v_Thist] = gradescent(@conformalcost,imax,aC,bC,etolC,0,...
    reshape(v',[],1),isedge,elsq_T);
  v_T = reshape(v_Thist(:,end),3,[])';
  [M_T,L_T] = lapbel(v_T,f_T);
  D_T = eigvf(L_T,M_T,numeig);

% import wavefront object file (3)
elseif target_case == 3
  filename = 'spot';
  fid = fopen(['../meshes/' filename '.obj'],'rt');
  [v_T,f_T] = readwfobj(fid);
  [M_T,L_T] = lapbel(v_T,f_T);
  D_T = eigvf(L_T,M_T,numeig);

elseif target_case == 4
  v_T = v - repmat(sphar(v),1,3).*vn;
  f_T = f;
  [M_T,L_T] = lapbel(v_T,f_T);
  D_T = eigvf(L_T,M_T,numeig);
end
%% initial conformal factors
s0 = exp(-zeros(numv,1));
%% MIEP2 via gradient descent
[J_hist,s] = gradescent(@eigencost,imax,aS,bS,etolS,0,...
  s0,M,L,D_T,numeig);
%% descent results
s_end = s(:,end);
D_endp = eigvf(L,diag(1./s_end)*M,numeig);
%% mfg. sol'n benchmark
if target_case == 2
  disp 'percent error'
  s_error = norm(s_T - s_end)./norm(s_T)*100
end
%% conformal fit
conf = sqrt(kron(1./s_end',1./s_end));
% elsq_end = elsq0.*conf;
elsq_end = elsq0.*conf(isedge); % linear indices
[Jc_hist,vhist] = gradescent(@conformalcost,imax,aC,bC,etolC,0,...
  reshape(v',[],1),isedge,elsq_end);
%% fit results
v_end = reshape(vhist(:,end),3,[])';
[M_end,L_end] = lapbel(v_end,f);
D_end = eigvf(L_end,M_end,numeig);
%% visualisation
close all;
% compare mesh
Mesh0 = TriRep(f,v); %triangulation(f,v);
Mesh_T = TriRep(f_T,v_T); %triangulation(f_T,v_T);
Mesh_end = TriRep(f,v_end); %triangulation(f,v_end);

figure(); set(gcf,'outerposition',[0, 0, 1024, 768]);
subplot(2,3,1); hold all; view(3); grid on; axis equal
trimesh(Mesh0);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('original mesh');

subplot(2,3,2); hold all; view(3); grid on; axis equal
trimesh(Mesh_T);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('target mesh');

subplot(2,3,3); hold all; view(3); grid on; axis equal
trimesh(Mesh_end);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('resultant mesh');

% compare spectra
subplot(2,3,4:6); hold all;
plot(D_0,'bx:');
if target_case == 2
  plot(D_Tp,'ks:');
end
plot(D_T,'k.-');
plot(D_endp,'rs:');
plot(D_end,'ro-');
if target_case == 2
  legend('\lambda_{initial}','\lambda_{Target pre}','\lambda_{Target}',...
    '\lambda_{result pre}','\lambda_{result}',...
    'location','best');
else
  legend('\lambda_{initial}','\lambda_{Target}',...
    '\lambda_{result pre}','\lambda_{result}',...
    'location','best');
end
xlabel('Number of eigenvalues'); ylabel('Eigenvalues of M^{-1}L');