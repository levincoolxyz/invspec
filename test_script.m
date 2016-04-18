clear;
%% some spherical harmonics
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
Y10 = @(v) v(:,3)./vnorm(v);
Y20 = @(v) (2*v(:,3).^2-v(:,1).^2-v(:,2).^2)./(vnorm(v)).^2;
Y32 = @(v) (v(:,1).^2-v(:,2).^2).*v(:,3)./(vnorm(v)).^3;
Y33 = @(v) (v(:,1).^2-3*v(:,2).^2).*v(:,1)./(vnorm(v)).^3;
Y43 = @(v) (7*v(:,3).^2-3*vnorm(v).^2).*v(:,1).*v(:,3)./(vnorm(v)).^4;
%% control parameters
imax = 3e3; % gradient descent maximum iterations
aC = .5; bC = .2; tC = 30; etolC = 5e-4; % Conformal descent control
aS = .5; bS = .4; tS = 150; etolS = 5e-4; % invSpec descent control
numeig = .6; % number of eigenvalues used, <1 => percent, <=0 => all
pert = .512; % scaling coefficient used to control target perturbation
rng(1432543); % rand seed
%% input case == 1; import face-vtx from *.obj file
% init_data.num = 1; 
% init_data.dat = 'sphere_small';
%% input case == 2; sphere of ssize # of vtx
init_data.num = 2; 
init_data.dat = '200';
%% input case == 3; import face-vtx from *.mat file
% init_data.num = 3; 
% init_data.dat = 'sphere500';
%% target case == 1; random conformal factor (on vtx) deformation
% target_data.num = 1;
%% target case == 2; prescribed perturbation (of sphere) along vtx normal
target_data.num = 2;
target_data.dat = @(v) abs(Y32(v));
%% target case == 3; import face-vtx from *.obj file
% target_data.num = 3;
% target_data.dat = 'spot';
for pert = linspace(.5,.7,3)
%% testing time
[v,f,v_end,v_T,f_T,J_hist,Jc_hist,...
  D_0,D_T,D_endp,D_end] = main(init_data,target_data,...
  imax,aC,bC,tC,etolC,aS,bS,tS,etolS,...
  numeig,pert);
%% visualing results
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
subplot(2,3,4:6); hold all; grid on;
plot((D_0 - D_T),'bx:');
% plot((D_Tp - D_T),'gs');
plot((D_endp - D_T),'k-','linewidth',2);
plot((D_end - D_T),'ro');
if exist('D_Tp','var')
  legend('\lambda_{initial} - \lambda_{target}',...
    '\lambda_{target embed} - \lambda_{target}',...
    '\lambda_{MIEP2} - \lambda_{target}',...
    '\lambda_{final embed} - \lambda_{target}',...
    'location','best');
else
  legend('\lambda_{initial} - \lambda_{target}',...
    '\lambda_{MIEP2} - \lambda_{target}',...
    '\lambda_{final embed} - \lambda_{target}',...
    'location','best');
end
% set(gca,'yscale','log');
xlabel('# of eigenvalues (#1 is of the highest frequency)'); 
ylabel('\propto eigenvalue magnitude');
title('Deviation from target Laplacian eigenvalues in magnitude');

ym = get(gca,'ylim');
text(floor(numeig/4.5),max(ym(2) + .18*diff(ym)),...
  num2str([J_hist(end) Jc_hist(end)],...
  ['Convergence Energies: J_{MIEP2} = %g     ',...
  'J_{embedding} = %g']));

%% store for record
if isa(target_data.dat,'function_handle')
  dumb = func2str(target_data.dat);
  dumb = dumb(5:end);
else
  dumb = target_data.dat;
end
endname = num2str([init_data.num, target_data.num, numeig, pert],...
  ['i%d_' init_data.dat '_t%d_' dumb '_e%gp%g']);
saveas(gcf,[endname '.png']);
save([endname '.mat'],'v','f','v_end','v_T','f_T',...
  'D_0','D_T','D_endp','D_end','J_hist','Jc_hist');
end
