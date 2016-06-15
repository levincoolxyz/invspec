clear;
load i2_300_t2_abs(Y33(v))_e0.1p0.5.mat
rng(1432543); % rand seed
numv = size(v,1);
[M,L] = lapbel(v,f);
[M_T,L_T] = lapbel(v_T,f_T);
D_Tf = eigvf(L_T,M_T,numv);

s = .5*rand(numv,1)+.75;
s = log(s);
D_Sf = eigvf(L,diag(1./s)*M,numv);
%%
% for h = [eps^.25 1e-6 sqrt(eps) eps]
for h = eps^.25
[~,GJ] = eigencost(s,M,L,D_T,30);
dJ = zeros(size(GJ));
% dJim = zeros(size(GJ));
for i = 1:numel(s)
  ds = zeros(size(s));
  ds(i) = h;
  Jph = eigencost(s+ds,M,L,D_T,30);
  Jmh = eigencost(s-ds,M,L,D_T,30);
  dJ(i) = (Jph-Jmh)/2/ds(i);
%   dJim(i) = imag(eigencost(s+complex(0,ds),M,L,D_T,20))/h;
end
%%
reldiff = (dJ-GJ)./dJ;
% reldiffim = (dJim-GJ)./dJ;
% mean(reldiff)
% std(reldiff)
%%
close all;
figure(); set(gcf,'outerposition',[0, 0, 1024, 768]);
subplot(2,2,1); hold all;
plot(reldiff,'kx--')
% plot(reldiffim,'ro--')
title('(\Delta J-\nabla J)/\Delta J');
subplot(2,2,2); hold all;
plot(GJ./dJ,'kx--')
title('\nabla J./\Delta J');
subplot(2,2,3:4); hold all;
plot(abs(dJ),'rx:')
plot(abs(GJ),'b-')
set(gca,'yscale','log')
legend('|\Delta J|','|\nabla J|');
saveas(gcf,num2str(h,'gradient verification h=%g.png'));
end