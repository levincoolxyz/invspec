clear; close all;
%% some spherical harmonics
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
Y10 = @(v) v(:,3)./vnorm(v);
Y20 = @(v) (2*v(:,3).^2-v(:,1).^2-v(:,2).^2)./(vnorm(v)).^2;
Y32 = @(v) (v(:,1).^2-v(:,2).^2).*v(:,3)./(vnorm(v)).^3;
Y33 = @(v) (v(:,1).^2-3*v(:,2).^2).*v(:,1)./(vnorm(v)).^3;
Y43 = @(v) (7*v(:,3).^2-3*vnorm(v).^2).*v(:,1).*v(:,3)./(vnorm(v)).^4;
%% control parameters
imax = 4e3; % descent maximum iteration
aS = .7; bS = .7; tS = 10; etolS = 1e-5; % MIEP2 descent control
aC = .5; bC = .8; tC = 10; etolC = 1e-4; % embedding descent control
maxL = 30; % max degree of SH basis considered
numeig = 16^2; % number of eigenvalues used, 0<x<=1 ratio, x<=0 full
numa = numeig; % just for the moment
% numa = 20^2; % number of free (nonzero) SH basis coefficients
pert = 0; % scaling coefficient used to control target perturbation
rng(1432543); % rand seed
method = 'BFGS'; % BFGS => fminunc, GD => in-house gradient descent
%% input case == 1; import face-vtx from *.obj file
% init_data.num = 1;
% init_data.dat = 'sphere_small';
%% input case == 2; sphere of ssize # of vtx
init_data.num = 2;
init_data.dat = '540';
% init_data.dat = '1000';
% init_data.dat = num2str(numeig);
%% input case == 3; import face-vtx from *.mat file [need v and f]
% init_data.num = 3;
% init_data.dat = '550';
%% input case == 4; use steady state cMCF of the target mesh
% init_data.num = 4;
% init_data.dat = 'mcf';
%% target case == 1; random vertex conformal factor deformations
% target_data.num = 1;
%% target case == 2; prescribed perturbation (of sphere) along vtx normal
% target_data.num = 2;
% target_data.dat = @(v) abs(Y32(v));
% target_data.dat = @(v) abs(Y43(v));
%% target case == 3; import face-vtx from *.obj file
target_data.num = 3;
% target_data.dat = 'bunny';
% target_data.dat = 'bunny1k';
% target_data.dat = 'blob1k';
target_data.dat = 'blob18k';
% target_data.dat = 'spot';
%% target case == 4; import face-vtx from *.mat file [need v_T and f_T]
% target_data.num = 4;
% target_data.dat = 'cow05';
%% target case == 5; prescribed spectra
% target_data.num = 5;
% target_data.name = 'bunnySH';
% load(num2str(sqrt(numeig)-1,'/home/ultimate/invspec/mcode/SH/bunnySHspecL=%d.mat'));
% target_data.a = a_pj;
% target_data.dat = eigvfSH(a_pj,numeig);
% target_data.name = 'sphere';
% D_s = [];
% for l = 0:20
%     D_s = [repmat(-l*(l+1),1,2*l+1), D_s];
% end
% target_data.dat = D_s';
%% testing time
% for pert = [.5:.1:.8 1 1.5 2]
% for reg = [1e-4 1e-3 1e-2 5e-2 .1 .5 1]
% for maxL = [10 20]
  if isa(target_data.dat,'function_handle')
    dumb = func2str(target_data.dat);
    dumb = dumb(5:end);
  elseif target_data.num == 5
    dumb = target_data.name;
  else
    dumb = target_data.dat;
  end
  if pert == 0
    endname = num2str([init_data.num, target_data.num, numa, numeig, maxL],...
      ['i%d_' init_data.dat '_t%d_' dumb '_a%ge%gL%g']);
  else
    endname = num2str([init_data.num, target_data.num, numa, numeig, maxL, pert],...
      ['i%d_' init_data.dat '_t%d_' dumb '_a%ge%gL%gp%g']);
  end
  %% main computation
  diary(['SH/' endname '.out']);
  [v,v_T,v_end,f,f_T,a_end,s_end,s_T,J_hist,Jc_hist,...
    D_0,D_T,D_endp,D_end,vmcf] = mainSH(init_data,target_data,...
    method,imax,aC,bC,tC,etolC,aS,bS,tS,etolS,...
    numeig,maxL,numa,pert);
  diary off;
  %% visualing results
  figh = visualizeSH(v,v_T,v_end,vmcf,f,f_T,s_end,s_T,...
    J_hist,Jc_hist,D_0,D_T,D_endp,D_end,0);
  %% record said results
  hgexport(figh,['SH/' endname '.png'],...
    hgexport('factorystyle'), 'Format', 'png'); 
  save(['SH/' endname '.mat'],'v','v_T','v_end','f','f_T','s_end','s_T',...
    'D_0','D_T','D_endp','D_end','J_hist','Jc_hist','a_end');
% end
