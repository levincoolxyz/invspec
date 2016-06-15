clear;
%% some spherical harmonics
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
Y10 = @(v) v(:,3)./vnorm(v);
Y20 = @(v) (2*v(:,3).^2-v(:,1).^2-v(:,2).^2)./(vnorm(v)).^2;
Y32 = @(v) (v(:,1).^2-v(:,2).^2).*v(:,3)./(vnorm(v)).^3;
Y33 = @(v) (v(:,1).^2-3*v(:,2).^2).*v(:,1)./(vnorm(v)).^3;
Y43 = @(v) (7*v(:,3).^2-3*vnorm(v).^2).*v(:,1).*v(:,3)./(vnorm(v)).^4;
%% control parameters
% imax = 1e4; % gradient descent maximum iterations
imax = 5e3; % gradient descent maximum iterations
aC = .5; bC = .8; tC = 10; etolC = 1e-7; % Conformal descent control
aS = .7; bS = .7; tS = 10; etolS = 1e-10; % invSpec descent control
% aC = .4; bC = .7; tC = 10; etolC = 1e-4; % Conformal descent control
% aS = .5; bS = .7; tS = 150; etolS = 1e-4; % invSpec descent control
numeig = 7; % number of eigenvalues used, 0<x<=1 percent, x<=0 full
pert = .5; % scaling coefficient used to control target perturbation
rng(1432543); % rand seed
%% input case == 1; import face-vtx from *.obj file
% init_data.num = 1; 
% init_data.dat = 'sphere_small';
%% input case == 2; sphere of ssize # of vtx
init_data.num = 2; 
init_data.dat = '300';
%% input case == 3; import face-vtx from *.mat file
% init_data.num = 3; 
% init_data.dat = '300';
%% target case == 1; random conformal factor (on vtx) deformation
% target_data.num = 1;
%% target case == 2; prescribed perturbation (of sphere) along vtx normal
% target_data.num = 2;
% target_data.dat = @(v) abs(Y32(v));
% target_data.dat = @(v) abs(Y33(v));
%% target case == 3; import face-vtx from *.obj file
target_data.num = 3;
% target_data.dat = 'bunny';
target_data.dat = 'bunny326';
% target_data.dat = 'spot487';
%% target case == 4; prescribed eigenvalue target
% target_data.num = 4;
% target_data.dat = 'D-T';
% D_s = [];
% for l = 0:20
%     D_s = [repmat(-l*(l+1),1,2*l+1), D_s];
% end
% target_data.D_T = D_s';
%% testing time
% for pert = [.5:.1:.8 1 1.5 2]
% for numeig = [.016 .019 .022 .025 .028 .032 .034 .037 .04 .044 .048 .064 .08 .1:.1:1]
  if isa(target_data.dat,'function_handle')
    dumb = func2str(target_data.dat);
    dumb = dumb(5:end);
  else
    dumb = target_data.dat;
  end
  endname = num2str([init_data.num, target_data.num, numeig, pert],...
    ['i%d_' init_data.dat '_t%d_' dumb '_e%gp%g']);
%  if target_data.num == 3
%    endname = [endname num2str(target_data.Nmcf,'_Nmcf%d')];
%  end
  diary([endname '.out']);
  %% main computation
  [v,v_T,v_end,f,f_T,s_end,s_T,J_hist,Jc_hist,...
    D_0,D_T,D_endp,D_end] = main(init_data,target_data,...
    imax,aC,bC,tC,etolC,aS,bS,tS,etolS,...
    numeig,pert);
  %% visualing results
  close all;
  figh = visualize(v,v_T,v_end,f,f_T,s_end,s_T,...
    J_hist,Jc_hist,D_0,D_T,D_endp,D_end);
  %% store for record
  hgexport(figh,[endname '.png'],...
    hgexport('factorystyle'), 'Format', 'png'); 
  save([endname '.mat'],'v','v_T','v_end','f','f_T','s_end','s_T',...
    'D_0','D_T','D_endp','D_end','J_hist','Jc_hist');
  diary off;
% end
