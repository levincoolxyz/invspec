clear 
% numeig = 602;
% filename = {'bunny','bunny327','bunny602','bunny1k','bunny2k','bunny4k','bunny7k'};

numeig = 1048;
filename = {'spot','spot487','spot1k'};

D = zeros(numeig,numel(filename));
Dmcf = zeros(numeig,numel(filename));
for i = 1:numel(filename)
  fid = fopen(['../meshes/' filename{i} '.obj'],'rt');
  [v_T,f_T] = readwfobj(fid);
  [M_T,L_T] = lapbel(v_T,f_T);
  load(['mcf/' filename{i} '.mat'])
  [M,L] = lapbel(v,f_T);
  if size(v_T,1) < numeig
    D(numeig-size(v_T,1)+1:end,i) = eigvf(L_T,M_T,size(v_T,1));
    Dmcf(numeig-size(v_T,1)+1:end,i) = eigvf(L,diag(1./s_T)*M,size(v_T,1));
  else
    D(:,i) = eigvf(L_T,M_T,numeig);
    Dmcf(:,i) = eigvf(L,diag(1./s_T)*M,numeig);
  end
end

%% plot results
D(D==0) = nan;
Dmcf(Dmcf==0) = nan;
close all;
rei = numeig:-1:1;
figure(); hold all; grid on;
for i = 1:numel(filename)
  plot(rei,-D(:,i));
end
for i = 1:numel(filename)
  plot(rei,-Dmcf(:,i),'--');
end
set(gca,'xlim',[1 numeig])
ylim([0 -min(min(D))]);
legend(filename,'location','northwest')
set(gcf,'outerposition',[0, 0, 1920, 1080]);

%% save results
% hgexport(gcf,'bunspec.png',hgexport('factorystyle'), 'Format', 'png'); 
% save('bunspec.mat','D','Dmcf');

hgexport(gcf,'cowspec.png',hgexport('factorystyle'), 'Format', 'png'); 
save('cowspec.mat','D','Dmcf');