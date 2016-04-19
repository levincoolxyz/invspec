close all;
figure(); hold all; grid on;
for pert = [0.4 .6 .7 .8 .9 1 1.5]
  numeig = ceil(.6*200);
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