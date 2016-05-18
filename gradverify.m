clear;
load i2_400_t3_bunny326_e20p0_Nmcf5000.mat
% load i2_200_t2_abs(Y33(v))_e0.1p5.mat
rng(1432543); % rand seed
numv = size(v,1);
[M,L] = lapbel(v,f);
[M_T,L_T] = lapbel(v_T,f_T);
%%
s = rand(numv,1);
[J,GJ] = eigencost(s,M,L,D_T,20);
dJ = zeros(size(GJ));
for i = 1:numel(s)
  ds = zeros(size(s));
  ds(i) = 1e-6;
  Ji = eigencost(s+ds,M,L,D_T,20);
  dJ(i) = (Ji-J)/ds(i);
end
%%
reldiff = (dJ-GJ)./dJ;
mean(reldiff)
std(reldiff)
%%
close all;
figure();
subplot(2,2,1); hold all;
plot(reldiff,'kx--')
title('(\Delta J-\nabla J)/\Delta J');
subplot(2,2,2); hold all;
plot(GJ./dJ,'kx--')
title('\nabla J./\Delta J');
subplot(2,2,3:4); hold all;
plot(dJ,'rx:')
plot(GJ,'b-')
set(gca,'yscale','log')
legend('\Delta J','\nabla J');