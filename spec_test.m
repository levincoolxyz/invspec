clear;
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
c = [.3 0 0];
%% theoretical spec
numeig = 1024;
N = 5;
D = nan(numeig,N+1);
Dinv = nan(numeig,N+1);

D_anal = [];
for l = 0:31
    D_anal = [repmat(-l*(l+1),1,2*l+1), D_anal];
end
D_anal = fliplr(D_anal);
%% do fine mesh
filename = 'sphere_large';
fid = fopen(['../meshes/' filename '.obj'],'rt');
[v,f] = readwfobj(fid);
scl0 = makeUnitArea(v,f);
[M,L] = lapbel(v,f);
D(:,1) = eigvf(L,M,numeig);
%% do inverted fine mesh
vinv = sphinv(v,c);
scl = makeUnitArea(vinv,f);
vinv = vinv - repmat(volCenter(vinv,f),size(v,1),1);
vinv = vinv*scl/scl0;
[Minv,Linv] = lapbel(vinv,f);
Dinv(:,1) = eigvf(Linv,Minv,numeig);
%% do coarse mesh with inverses and plot spectrum
for i = 2:N+1
  ssize = 2^(i-2)*100;
  if exist([num2str(ssize,'%g') '.mat'],'file')
    load([num2str(ssize,'%g') '.mat']);
  else
    v = ParticleSampleSphere('Vo',RandSampleSphere(ssize));
    f = fliplr(convhulln(v));
    save([num2str(ssize,'%g') '.mat'],'f','v');
  end
  [M,L] = lapbel(v,f);
  scl0 = makeUnitArea(v,f);
  vinv = sphinv(v,c);
  scl = makeUnitArea(vinv,f);
  vinv = vinv - repmat(volCenter(vinv,f),size(v,1),1);
  vinv = vinv*scl/scl0;
  [Minv,Linv] = lapbel(vinv,f);
  if ssize<numeig
    D(numeig-ssize+1:numeig,i) = eigvf(L,M,numeig);
    Dinv(numeig-ssize+1:numeig,i) = eigvf(Linv,Minv,numeig);
  else
    D(:,i) = eigvf(L,M,numeig);
    Dinv(:,i) = eigvf(Linv,Minv,numeig);
  end
end
D = flipud(D);
Dinv = flipud(Dinv);
%% setting up the graph
close all;
figure(); hold all; grid on;
for i = 1:N+1
  plot(abs((D(:,i)-D_anal')./D_anal'),'-');
  plot(abs((Dinv(:,i)-D_anal')./D_anal'),'--');
end
ym = get(gca,'ylim');
set(gca,'xlim',[2 numeig],'ylim',[0 ym(2)],'xscale','log','yscale','log');
legend('40962','40962inv','100','100inv','200','200inv','400','400inv',...
  '800','800inv','1600','1600inv','location','southeast');
title('Relative eigenvalue differences');
xlabel('# of eigenvalues');
saveas(gcf,'sphere_spec_test.png');

%% visualize spherical inversion on sphere
scl0 = makeUnitArea(v,f);
vinv = sphinv(v,c);
% vinv = sphinv(v,[.5 0 0]);
% vinv = sphinv(v,[.9 0 0]);
scl = makeUnitArea(vinv,f);
vinv = vinv - repmat(volCenter(vinv,f),size(v,1),1);
vinv = vinv*scl/scl0;

% figure(); hold all; view(3); grid on; axis equal
% trimesh(f,v(:,1),v(:,2),v(:,3));
figure(); hold all; view(3); grid on; axis equal
trimesh(f,vinv(:,1),vinv(:,2),vinv(:,3));