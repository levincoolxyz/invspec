clear;
load '/home/ultimate/invspec/mcode/SH/bunnySHspecL=30.mat'
rng(1432543); % rand seed
numa = 31^2;
D_T = eigvfSH(a_pj,numa);

% a = [7 2 .2 -.4 2 0 1 0 -1 -1 0 .5 zeros(1,numa-12)]' + .5*rand(numa,1)./(1:numa)';
a = [2*sqrt(pi);zeros(numa-1,1)] + .5*rand(numa,1)./(1:numa)';
numeig = numa;
%%
% for h = [eps^.25 1e-6 sqrt(eps) eps]
for h = 1e-6
dJ = zeros(numa,1);
% dJim = zeros(size(GJ));
for i = 1:numel(a)
  da = zeros(size(a));
  da(i) = h;
  Jph = eigencostSH(a+da,D_T,numeig);
  Jmh = eigencostSH(a-da,D_T,numeig);
  dJ(i) = (Jph-Jmh)/2/da(i);
end
%%
[~,GJ] = eigencostSH(a,D_T,numeig);
reldiff = (dJ-GJ)./dJ;
% reldiffim = (dJim-GJ)./dJ;
% mean(reldiff)
% std(reldiff)
%%
% close all;
figure(); set(gcf,'outerposition',[0, 0, 1024, 768]);
subplot(2,2,1); hold all;
plot(reldiff,'kx--')
% plot(reldiffim,'ro--')
title('(\Delta J-\nabla J)/\Delta J');
subplot(2,2,2); hold all;
plot(GJ./dJ,'kx--')
title('\nabla J./\Delta J');
subplot(2,2,3:4); hold all;
% plot((dJ),'rx:')
% plot((GJ),'b-')
plot(abs(dJ),'rx:')
plot(abs(GJ),'b-')
set(gca,'yscale','log')
legend('|\Delta J|','|\nabla J|');
% saveas(gcf,num2str(h,'SH/gradient verification h=%g.png'));
end