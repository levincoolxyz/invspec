close all; clear;
%% initialization
input_case = 2; % 1 - octa; 2 - icosa; 3 - obj file
target_case = 2; % 1 - octa; 2 - rand conf deform; 3 - obj file
alpha = .5; beta = .8; imax = 100; % gradient descent control param
numeig = 0; % number of eigenvalues used, 0 means full input
rng(1432543); % rand seed
def = .6; % scaling coefficient used for target case #2
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
  temp = TriQuad({f v});
  f = temp{1};
  v = temp{2};
  
% import wavefront object file (3)
elseif input_case == 3
  filename = 'sphere_small';
  fid = fopen(['../meshes/' filename '.obj'],'rt');
  [v,f] = readwfobj(fid);
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
%% compute initial edge lengths squared
elsq0 = zeros(numv);
for i = 1:numv
  for j = (i+1):numv % skipping the symmetric lower triangular part
    elsq0(i,j) = sum((v(i,:)-v(j,:)).^2);
  end
end
elsq0 = elsq0.*isedge;
%% compute laplacian
[M,L] = lapbel(v,f);
% D_0 = sort(eig(inv(M)*L));
% D_0 = sort(eigs(L,M,numeig,0));
D_0 = eigvf(L,M,numeig);
%% target spectrum
% skewed octahedron
if target_case == 1
  v_T = [0 0 .5; 0 1 0; 2 0 0 ;-2 0 0; 0 -3 0; 0 0 -.5];
  f_T = fliplr(convhulln(v));
  [M_T,L_T] = lapbel(v_T,f_T);
%   D_T = sort(eig(inv(M_T)*L_T));
  D_T = eigvf(L_T,M_T,numeig);

% random small prescribed conformal deformation
elseif target_case == 2
  s_T = exp(-rand(numv,1)*def);
%   D_Tp = sort(eig(diag(s_T)*inv(M)*L));
  D_Tp = eigvf(L,diag(1./s_T)*M,numeig);
  f_T = f;
  conf_T = sqrt(kron(1./s_T',1./s_T));
  elsq_T = elsq0.*conf_T;
  [Jc_Thist,vhist_T] = gradescent(@conformalcost,imax,alpha,beta,...
    reshape(v',[],1),isedge,elsq_T);
  Jc_T = Jc_Thist(end)
  v_T = reshape(vhist_T(:,end),3,[])';
  [M_T,L_T] = lapbel(v_T,f_T);
%   D_T = sort(eig(inv(M_T)*L_T));
  D_T = eigvf(L_T,M_T,numeig);

% import wavefront object file (3)
elseif input_case == 3
  filename = 'blub';
  fid = fopen(['../meshes/' filename '.obj'],'rt');
  [v_T,f_T] = readwfobj(fid);
  [M_T,L_T] = lapbel(v_T,f_T);
%   D_T = sort(eig(inv(M_T)*L_T));
  D_T = eigvf(L_T,M_T,numeig);
  D_T = D_T((end-numeig):end);
end
%% initial conformal factors
s0 = exp(-zeros(numv,1));
%% LiIEP via gradient descent
[J,s] = gradescent(@eigencost,imax,alpha,beta,s0,M,L,D_T,numeig);
%% descent convergence
% figure(); hold all; grid on;
% plot(J,'k.:')
% set(gca,'xscale','log','yscale','log'); axis square
% % linearregress(log(1:imax),log(J),1);
%% descent results
J_end = J(end)
s_end = s(:,end);
% D_endp = sort(eig(diag(s_end)*inv(M)*L));
D_endp = eigvf(L,diag(1./s_end)*M,numeig);
%% mfg. sol'n benchmark
if target_case == 2
  disp 'percent error'
  s_error = norm(s_T - s_end)./norm(s_T)*100
end
%% conformal fit
conf = sqrt(kron(1./s_end',1./s_end));
elsq_end = elsq0.*conf;
[Jc,vhist] = gradescent(@conformalcost,imax,alpha,beta,...
  reshape(v',[],1),isedge,elsq_end);
%% fit convergence
% figure(); hold all; grid on;
% plot(Jc,'k.:')
% set(gca,'xscale','log','yscale','log'); axis square
% % linearregress(log(1:imax),log(J),1);
%% fit results
Jc_end = Jc(end)
v_end = reshape(vhist(:,end),3,[])';
[M_end,L_end] = lapbel(v_end,f);
% D_end = sort(eig(inv(M_end)*L_end));
D_end = eigvf(L_end,M_end,numeig);
%% visualisation
% difference in mesh
Mesh0 = TriRep(f,v); %triangulation(f,v);
Mesh_T = TriRep(f_T,v_T); %triangulation(f_T,v_T);
Mesh_end = TriRep(f,v_end); %triangulation(f,v_end);

figure(); set(gcf,'outerposition',[0, 0, 1024, 768]);
subplot(2,3,1); hold all; view(3); grid on;
trimesh(Mesh0);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('original mesh');

subplot(2,3,2); hold all; view(3); grid on;
trimesh(Mesh_T);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('target mesh');

subplot(2,3,3); hold all; view(3); grid on;
trimesh(Mesh_end);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('resultant mesh');

% difference in spectra
subplot(2,3,4:6); hold all;
plot(D_0,'bx:');
plot(D_Tp,'ks:');
plot(D_T,'k.-');
plot(D_endp,'rs:');
plot(D_end,'ro-');
legend('\lambda_{initial}','\lambda_{Target pre}','\lambda_{Target}',...
  '\lambda_{result pre}','\lambda_{result}',...
  'location','best');
xlabel('Number of eigenvalues'); ylabel('Eigenvalues of M^{-1}L');