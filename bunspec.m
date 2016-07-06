numeig = 602;
filename = {'bunny','bunny327','bunny602','bunny1k','bunny2k','bunny4k','bunny7k'};
for i = 1:7
  fid = fopen(['../meshes/' filename{i} '.obj'],'rt');
  [v,f] = readwfobj(fid);
  [M,L] = lapbel(v,f);
  if size(v,1) <= numeig
    D(numeig-size(v,1)+1:end,i) = eigvf(L,M,size(v,1));
  else
    D(:,i) = eigvf(L,M,numeig);
  end
end
%%
D(D==0) = nan;
close all;
rei = numeig:-1:1;
figure(); hold all; grid on;
for i = 1:7
  plot(rei,-D(:,i));
end
set(gca,'xlim',[1 numeig])
ylim([0 -min(min(D))]);
legend(filename,'location','northwest')
set(gcf,'outerposition',[0, 0, 1920, 1080]);
hgexport(gcf,'bunspec.png',...
  hgexport('factorystyle'), 'Format', 'png'); 