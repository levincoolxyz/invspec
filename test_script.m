clear;
%% some spherical harmonics
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
Y10 = @(v) v(:,3)./vnorm(v);
Y20 = @(v) (2*v(:,3).^2-v(:,1).^2-v(:,2).^2)./(vnorm(v)).^2;
Y32 = @(v) (v(:,1).^2-v(:,2).^2).*v(:,3)./(vnorm(v)).^3;
Y33 = @(v) (v(:,1).^2-3*v(:,2).^2).*v(:,1)./(vnorm(v)).^3;
Y43 = @(v) (7*v(:,3).^2-3*vnorm(v).^2).*v(:,1).*v(:,3)./(vnorm(v)).^4;
%% control parameters
% imax = 3e3; % gradient descent maximum iterations
% aC = .5; bC = .2; tC = 30; etolC = 5e-4; % Conformal descent control
% aS = .5; bS = .4; tS = 150; etolS = 5e-4; % invSpec descent control
imax = 1e3; % gradient descent maximum iterations
aC = .5; bC = .8; tC = 30; etolC = 5e-4; % Conformal descent control
aS = .5; bS = .8; tS = 150; etolS = 5e-4; % invSpec descent control
numeig = .6; % number of eigenvalues used, <=1 => percent, <=0 => all
pert = .512; % scaling coefficient used to control target perturbation
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
target_data.num = 2;
target_data.dat = @(v) abs(Y33(v));
%% target case == 3; import face-vtx from *.obj file
% target_data.num = 3;
% target_data.dat = 'spot';
%% testing time
% for pert = .3:.1:2
% for numeig = .1:.1:1
for numeig = .6
  [v,f,v_end,v_T,f_T,J_hist,Jc_hist,...
    D_0,D_T,D_endp,D_end] = main(init_data,target_data,...
    imax,aC,bC,tC,etolC,aS,bS,tS,etolS,...
    numeig,pert);
  %% visualing results
  close all;
  figh = visualize(v,f,v_end,v_T,f_T,...
    J_hist,Jc_hist,D_0,D_T,D_endp,D_end);
  %% store for record
  if isa(target_data.dat,'function_handle')
    dumb = func2str(target_data.dat);
    dumb = dumb(5:end);
  else
    dumb = target_data.dat;
  end
  endname = num2str([init_data.num, target_data.num, numeig, pert],...
    ['i%d_' init_data.dat '_t%d_' dumb '_e%gp%g']);
  saveas(figh,[endname '.png']);
  save([endname '.mat'],'v','f','v_end','v_T','f_T',...
    'D_0','D_T','D_endp','D_end','J_hist','Jc_hist');
end
