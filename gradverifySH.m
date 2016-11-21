clear;
load 'SH/forward/blob18kSHspecL=30.mat'
rng(1432543); % rand seed
% numa = 7^2;
numa = 10;
maxL = 7;
D_T = eigvfSH([a_pj(1:numa);zeros((maxL+1)^2-numa,1)],numa,maxL);

% a = [7 2 .2 -.4 2 0 1 0 -1 -1 0 .5 zeros(1,numa-12)]' + .5*rand(numa,1)./(1:numa)';
% a = [2*sqrt(pi);zeros(numa-1,1)] + .5*rand(numa,1)./(1:numa)';
a = [a_pj(1:numa) + .5*rand(numa,1)./(1:numa)';zeros((maxL+1)^2-numa,1)];
numeig = maxL^2;
%%
% for h = [eps^.25 1e-6 sqrt(eps) eps]
for h = 1e-6
dJ = zeros(numa,1);
for j = 1:numa
  da = zeros(size(a));
  da(j) = h;
  Jph = eigencostSH(a+da,D_T,numeig,maxL);
  Jmh = eigencostSH(a-da,D_T,numeig,maxL);
  dJ(j) = (Jph-Jmh)/2/h;
end
%%
[~,GJ] = eigencostSH(a,D_T,numeig,maxL);
GJ = GJ(1:numa);
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