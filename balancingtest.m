clear;
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
c0 = [.3 0 0];
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
vinv = sphinv(v,c0);
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
  vinv = sphinv(v,c0);
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

%% visualize double spherical inversion on sphere
scl0 = makeUnitArea(v,f);
vinv = sphinv(v,c0);
% vinv = sphinv(v,[.5 0 0]);
% vinv = sphinv(v,[.9 0 0]);
scl = makeUnitArea(vinv,f);
vinv = vinv - repmat(volCenter(vinv,f),size(v,1),1);
vinv = vinv*scl/scl0;

% figure(); hold all; view(3); grid on; axis equal
% trimesh(f,v(:,1),v(:,2),v(:,3));
figure(); hold all; view(3); grid on; axis equal
trimesh(f,vinv(:,1),vinv(:,2),vinv(:,3));
%% double spherical inversion test
clear
load mcf/blob1k.mat
numv = size(v_T,1);
c0 = [.1 .1 .1];
scl = makeUnitArea(v_T,f_T);
v_T = v_T - repmat(volCenter(v_T,f_T),numv,1);
v_T = v_T*scl*sqrt(4*pi);
vinv = sphinv(v_T,c0);
% vinv = sphinv(vinv,c0);
scl = makeUnitArea(vinv,f_T);
vinv = vinv - repmat(volCenter(vinv,f_T),numv,1);
vinv = vinv*scl*sqrt(4*pi);
vinv2 = sphinv(vinv,-c0);
scl = makeUnitArea(vinv2,f_T);
vinv2 = vinv2 - repmat(volCenter(vinv2,f_T),numv,1);
vinv2 = vinv2*scl*sqrt(4*pi);

close all; figure();
set(gcf,'outerposition',[0, 0, 1920, 540]);
 
subplot(1,3,1); hold all; view(3); grid on; axis equal
trimesh(f_T,v_T(:,1),v_T(:,2),v_T(:,3));
subplot(1,3,2); hold all; view(3); grid on; axis equal
trimesh(f_T,vinv(:,1),vinv(:,2),vinv(:,3));
subplot(1,3,3); hold all; view(3); grid on; axis equal
trimesh(f_T,vinv2(:,1),vinv2(:,2),vinv2(:,3));
%% mobius balancing test
clear
load mcf/blob1k.mat
% fid = fopen('../meshes/bunny1043.obj','rt');
% [v_T,f_T] = readwfobj(fid);

numv = size(v_T,1);
c0 = [.5 0 .1];
scl = makeUnitArea(v_T,f_T);
v_T = v_T - repmat(volCenter(v_T,f_T),numv,1);
v_T = v_T*scl*sqrt(4*pi);
vinv = sphinv(v_T,c0);
scl = makeUnitArea(vinv,f_T);
vinv = vinv - repmat(volCenter(vinv,f_T),numv,1);
vinv = vinv*scl*sqrt(4*pi);
% vinv = sphinv(vinv,[0 .3 .2]);
% scl = makeUnitArea(vinv,f_T);
% vinv = vinv - repmat(volCenter(vinv,f_T),numv,1);
% vinv = vinv*scl*sqrt(4*pi);

[vba1,c1] = mobiusbalancing(v_T,v_T,f_T);
scl = makeUnitArea(vba1,f_T);
vba1 = vba1 - repmat(volCenter(vba1,f_T),numv,1);
vba1 = vba1*scl*sqrt(4*pi);
[vba2,c2] = mobiusbalancing(vinv,v_T,f_T);
scl = makeUnitArea(vba2,f_T);
vba2 = vba2 - repmat(volCenter(vba2,f_T),numv,1);
vba2 = vba2*scl*sqrt(4*pi);

% dcm = pa(v_T,f_T); % princomp(v_T);
% % dcm = [dcm(:,2) dcm(:,3) dcm(:,1)];
% dcm = dcm*det(dcm); % ensure orientation preserving
% v_T = dcmrot(v_T,dcm);
% 
% dcm = pa(vba2,f_T); % princomp(v_T);
% dcm = dcm*det(dcm); % ensure orientation preserving
% vba2 = dcmrot(vba2,dcm);

close all; figure();
set(gcf,'outerposition',[0, 0, 1920, 1080]);
 
subplot(2,2,1); hold all; view(3); grid on; axis equal
title('original mesh')
xlabel('x'); ylabel('y'); zlabel('z');
trimesh(f_T,v_T(:,1),v_T(:,2),v_T(:,3));
text(4,1,-2,num2str([0 0 0],'%g\n'));

subplot(2,2,3); hold all; view(3); grid on; axis equal
title('Moebius-balanced mesh')
xlabel('x'); ylabel('y'); zlabel('z');
trimesh(f_T,vba1(:,1),vba1(:,2),vba1(:,3));
text(4,1,-2,num2str(c1,'%g\n'));

subplot(2,2,2); hold all; view(3); grid on; axis equal
title('apply double inversion')
xlabel('x'); ylabel('y'); zlabel('z');
trimesh(f_T,vinv(:,1),vinv(:,2),vinv(:,3));
text(-4,1,0.5,num2str(c0,'%g\n'));

subplot(2,2,4); hold all; view(3); grid on; axis equal
title('Moebius-balanced inverted mesh')
xlabel('x'); ylabel('y'); zlabel('z');
trimesh(f_T,vba2(:,1),vba2(:,2),vba2(:,3));
text(-4,1,0.5,num2str(c2,'%g\n'));