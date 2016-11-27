clear; close all; 
numalist = [7 10 16 17 20].^2;
lg = {};
figure(); hold all; grid on;
set(gcf,'outerposition',[0, 0, 1920, 1080]);
for numa = numalist
  endname = num2str(numa,'i2_540_t3_blob18k_a%de256L30');
  load(['SH/' endname '.mat']);
  a = nonzeros(a_end);
  L = 0:sqrt(numa)-1;
  amean = zeros(size(L));
  aup = zeros(size(L));
  adown = zeros(size(L));
  for l = L
    al = (abs(a((1:2*l+1)+l^2)));
    amean(l+1) = mean(al);
    aup(l+1) = max(al)-mean(al);
    adown(l+1) = mean(al)-min(al);
  end
  errorbar(L,log10(amean),.434*(adown./amean),.434*(aup./amean),...
    ':x','linewidth',2,'markersize',15)
  lg = [lg;num2str(numa)];
end

%% regularized
numalist = [18 20 25].^2;
for numa = numalist
  endname = num2str(numa,'i2_540_t3_blob18k_a%de256L30r*');
  endname = dir(['SH/' endname '.mat']);
  [~,endname] = fileparts(endname.name);
  load(['SH/' endname '.mat']);
%   a = nonzeros(a_end);
  a = a_end(1:numa);
  L = 0:sqrt(numa)-1;
  amean = zeros(size(L));
  aup = zeros(size(L));
  adown = zeros(size(L));
  for l = L
    al = (abs(a((1:2*l+1)+l^2)));
    amean(l+1) = mean(al);
    aup(l+1) = max(al)-mean(al);
    adown(l+1) = mean(al)-min(al);
  end
  errorbar(L,log10(amean),.434*(adown./amean),.434*(aup./amean),...
    '--o','linewidth',2,'markersize',12)
  lg = [lg;num2str(numa,'%dreg')];
end
%% tags and labels
plot(L,log10(1./L),'k-');
plot(L,log10(1./L.^2),'k-');
% plot(L,log10(exp(-L)),'k--');
legend([lg;'L^{-1}';'L^{-2}']);%;'exp(-L)']);
set(gca,'xlim',[0 max(L)]);
ylabel('Coefficient Magnitude $log(|a|)$','interpreter','latex');
xlabel('SH Degree $L$','interpreter','latex');
title('Spherical Harmonic Coefficients vs. its Degree');
hgexport(gcf,['SH/amag.png'],...
  hgexport('factorystyle'), 'Format', 'png'); 