close all;
figure(); hold all; grid on;
numeig = ceil(.6*200);
for pert = [0.4 .6 .7 .8 .9 1 1.5]
  load(num2str(pert,'i2_200_t2_abs(Y32(v))_e0.6p%g.mat'));
  [M,L] = lapbel(v,f);
  [M_T,L_T] = lapbel(v_T,f_T);
  [M_end,L_end] = lapbel(v_end,f);
  D_T = eigvf(L_T,M_T,numeig);
  %%
%   figure(); hold all; grid on;
%   ylim([0 .2]);
%   plot(diag(M))
%   plot(diag(M_T))
%   plot(diag(M_end))
  plot(D_T);
  pause(1);
end

%%
close all;
figure(); hold all; grid on;
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
Y32 = @(v) (v(:,1).^2-v(:,2).^2).*v(:,3)./(vnorm(v)).^3;
target_data.dat = @(v) abs(Y32(v));
numeig = 100;
pert = .6; % scaling coefficient used to control target perturbation
for ssize = 100:100:900
  v = ParticleSampleSphere('Vo',RandSampleSphere(ssize));
  f = fliplr(convhulln(v));
  [M,L] = lapbel(v,f);
  Hn = .5*[inv(M)*L*v(:,1) inv(M)*L*v(:,2) inv(M)*L*v(:,3)];
  H = vnorm(Hn);
  vn = Hn./repmat(H,1,3);
  v_T = v - repmat(target_data.dat(v),1,3).*vn*pert;
  f_T = f;
  [M_T,L_T] = lapbel(v_T,f_T);
  D_T = eigvf(L_T,M_T,numeig);
  plot(D_T);
  pause(1);
end

xlabel('# of eigenvalues'); ylabel('Eigenvalue of M^{-1}L');
title('First 100 eigenvalues of |Y32| perturbed sphere with different mesh refinement');
saveas(gcf,'Y32_spectrum.png');
save('Y32_spectrum_900.mat','D_T');
%% compare sphere spectrum
close all;
numeig = 100;
ssize = 900;
v = ParticleSampleSphere('Vo',RandSampleSphere(ssize));
f = fliplr(convhulln(v));
[M,L] = lapbel(v,f);

D_s = [];
for l = 0:9
    D_s = [repmat(-l*(l+1),1,2*l+1), D_s];
end
D = eigvf(L,M,numeig);
figure(); hold all; grid on;
plot(D_s)
plot(D)
legend('900 vtx sphere','theoretical','location','best');
xlabel('# of eigenvalues'); ylabel('Eigenvalue of M^{-1}L');
title('First 100 eigenvalues of sphere');
saveas(gcf,'sphere_spectrum_900.png');
%%
close all;
numeig = 100;
ssize = 100;
v = ParticleSampleSphere('Vo',RandSampleSphere(ssize));
f = fliplr(convhulln(v));
[M,L] = lapbel(v,f);

D_s = [];
for l = 0:9
    D_s = [repmat(-l*(l+1),1,2*l+1), D_s];
end
D = eigvf(L,M,numeig);
figure(); hold all; grid on;
plot(D_s)
plot(D)
legend('100 vtx sphere','theoretical','location','best');
xlabel('# of eigenvalues'); ylabel('Eigenvalue of M^{-1}L');
title('First 100 eigenvalues of sphere');
saveas(gcf,'sphere_spectrum_100.png');