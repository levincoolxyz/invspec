%% test on diff resolution mesh
close all;
target = 'blob';
res = {'1k','4k','18k','96k'};
maxL = 20;
numeig = (maxL+1)^2;

figure(); hold all; grid on;
for i = 1:numel(res)
load(num2str(maxL,['SH/forward/' target res{i} 'SHspecL=%d.mat']));
plot(numeig:-1:1,abs((-D_pj+D_T)./D_T),'--','linewidth',2);
end
set(gca,'yscale','log');
xlabel('# of Eigenvalues');
ylabel('Relative Error');
title(target);
legend(res,'location','best');
saveas(gcf,num2str(maxL,['SH/forward/' target 'SHspecLogDiffL=%d.png']));

%% test on diff maxL
close all;
% target = 'blob1k';
% Lrange = [1:3 5 8 10 15 20 21 22];
target = 'blob18k';
% target = 'blob96k';
Lrange = [1:3 5 8 10 15 20 25 30];

figure(); hold all; grid on;
for i = 1:numel(Lrange)
maxL = Lrange(i);
numeig = (maxL+1)^2;
load(num2str(maxL,['SH/forward/' target 'SHspecL=%d.mat']));
plot(numeig:-1:1,abs((-D_pj(end-numeig+1:end)+D_T)./D_T));
% plot(numeig:-1:1,abs((-D_ls(end-numeig+1:end)+D_T)./D_T));
end
set(gca,'yscale','log');
xlabel('# of Eigenvalues');
ylabel('Relative Error');
title(target);
legend(num2str(Lrange'),'location','best');
saveas(gcf,num2str(maxL,['SH/forward/' target 'SHspecLogDiffL=%d.png']));