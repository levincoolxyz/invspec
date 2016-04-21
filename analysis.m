
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
%% control parameters
imax = 3e3; % gradient descent maximum iterations
aC = .5; bC = .2; tC = 30; etolC = 5e-4; % Conformal descent control
aS = .5; bS = .4; tS = 150; etolS = 5e-4; % invSpec descent control
numeig = .6; % number of eigenvalues used, <=1 => percent, <=0 => all
pert = .512; % scaling coefficient used to control target perturbation
rng(1432543); % rand seed

%% test perturbation effect
% close all;
% figure(); hold all; grid on;
% numeig = ceil(.6*200);
% for pert = [0.4 .6 .7 .8 .9 1 1.5]
%   load(num2str(pert,'i2_200_t2_abs(Y32(v))_e0.6p%g.mat'));
%   [M,L] = lapbel(v,f);
%   [M_T,L_T] = lapbel(v_T,f_T);
%   [M_end,L_end] = lapbel(v_end,f);
%   D_T = eigvf(L_T,M_T,numeig);

% %   figure(); hold all; grid on;
% %   ylim([0 .2]);
% %   plot(diag(M))
% %   plot(diag(M_T))
% %   plot(diag(M_end))
%   plot(D_T);
%   pause(1);
% end

%% test mesh refinement laplacian eig convergence
close all;
figure(); hold all; grid on;
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
Y32 = @(v) (v(:,1).^2-v(:,2).^2).*v(:,3)./(vnorm(v)).^3;
target_data.dat = @(v) abs(Y32(v));
numeig = 100;
pert = .6; % scaling coefficient used to control target perturbation
for ssize = 1e2:1e2:1e3
  v = ParticleSampleSphere('Vo',RandSampleSphere(ssize));
  f = fliplr(convhulln(v));
  [M,L] = lapbel(v,f);
  D(:,ssize/100) = eigvf(L,M,numeig);
  Hn = .5*[inv(M)*L*v(:,1) inv(M)*L*v(:,2) inv(M)*L*v(:,3)];
  H = vnorm(Hn);
  vn = Hn./repmat(H,1,3);
  v_T = v - repmat(target_data.dat(v),1,3).*vn*pert;
  f_T = f;
  [M_T,L_T] = lapbel(v_T,f_T);
  D_T = eigvf(L_T,M_T,numeig);
  plot(D_T);
%   pause(1);
end

legend('100','200','300','400','500','600','700','800','900','1E3',...
  'location','best');
xlabel('# of eigenvalues'); ylabel('Eigenvalue of M^{-1}L');
title(['First 100 eigenvalues of |Y32| perturbed sphere' ...
  'with different mesh refinement']);
saveas(gcf,'Y32_spectrum.png');

%% compare mesh refinement on sphere spectrum
close all;
numeig = 1000;
D_s = [];
for l = 0:29
    D_s = [repmat(-l*(l+1),1,2*l+1), D_s];
end
filename = 'sphere_large';
fid = fopen(['../meshes/' filename '.obj'],'rt');
[v,f] = readwfobj(fid);
[M,L] = lapbel(v,f);
D = eigvf(L,M,numeig);
  
figure(); hold all; grid on;
plot(D)
plot(D_s)
legend('100','200','300','400','500','600','700','800','900','1E3',...
  'theoretical','location','best');
xlabel('# of eigenvalues'); ylabel('Eigenvalue of M^{-1}L');
title('First 100 eigenvalues of sphere');
saveas(gcf,'sphere_spectrum_fine.png');

%% compare sphere spectrum on fine mesh
close all;
D_s = [];
for l = 0:9
    D_s = [repmat(-l*(l+1),1,2*l+1), D_s];
end
figure(); hold all; grid on;
for i = 1:size(D,2)
  plot(D(:,i))
end
plot(D_s)
legend('sphere_large',...
  'theoretical','location','best');
xlabel('# of eigenvalues'); ylabel('Eigenvalue of M^{-1}L');
title('First 1000 eigenvalues of sphere');
saveas(gcf,'sphere_spectrum.png');

%% cMCF test
numeig = ceil(.6*200);
imax = 3e3;
% for pert = [0.4 .6 .7 .8 .9 1 1.5]
for pert = 1
  load(num2str(pert,'i2_200_t2_abs(Y32(v))_e0.6p%g.mat'));
  [M,L] = lapbel(v,f);
  [M_T,L_T] = lapbel(v_T,f_T);
  %% conformal Mean Curvature Flow for conformal factors
  s_T = meancurvflow(v_T,f_T,1000,'c',L_T,M_T);
  % D_Tp = eigvf(L,diag(1./s_T)*M,numeig);
  %% can I flow it back?
  conf_T = sqrt(kron(1./s_T',1./s_T));
  elsq_T = elsq0.*conf_T(isedge); % linear indices
  [~,v_Thist] = gradescent(@conformalcost,imax,aC,bC,tC,etolC,0,...
    reshape(v',[],1),isedge,elsq_T);
  v_c = reshape(v_Thist(:,end),3,[])';
  norm(vnorm(v_T - v_c))
  Mesh_T = TriRep(f_T,v_T); %triangulation(f_T,v_T);
  Mesh_c = TriRep(f,v_c); %triangulation(f,v_end);
end