clear;
%% load target and mcf data
% target = 'spot';
% target = 'spot10k';
target = 'bunny';
% target = 'spotS10';
% target = 'sphere_large';

% discriptor = '';
discriptor = 'inv';
% discriptor = 'log';

if exist(['mcf/' target '.mat'],'file')
  load(['mcf/' target '.mat'])
else
  fid = fopen(['../meshes/' target '.obj'],'rt');
  if fid == -1 % try .mat formats
    load(['../meshes/' target '.mat']);
  else
    [v_T,f_T] = readwfobj(fid);
  end
  [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');
  save(['mcf/' target '.mat'],'v_T','f_T','s_T','v');
end

% pert = @(x) .5*(2*v(:,3).^2-v(:,1).^2-v(:,2).^2)./(vnorm(v)).^2;
% 
% [M,L] = lapbel(v_T,f_T);
% invM = diag(1./diag(M));
% Hn = .5*[invM*L*v(:,1) invM*L*v(:,2) invM*L*v(:,3)];
% H = vnorm(Hn);
% vn = Hn./repmat(H,1,3);
% v_T = v_T - repmat(pert(v),1,3).*vn;
% [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');

if strcmp('inv',discriptor)
elseif strcmp('log',discriptor)
  s_T = log(1./s_T);
else
  s_T = 1./s_T; % re-invert / restore conf. factors
end

[vham,fham] = xyz2hammer(v);
%%
apjlist = [];
alslist = [];
Lrange = 10;%[1:3 5 8 10 15 20 25 30];
for maxL = Lrange
%% get spherical harmonics on vertices
[Y_v,LM] = sphericalHarmonicBase(v,maxL);
numSH = size(LM,1);

%% calculate conf. factors in SH basis
x = s_T;

a_ls = (Y_v'*Y_v)\(Y_v'*x);
s_ls = Y_v*a_ls;

M = lapbel(v,f_T);
delta = Y_v'*M*Y_v;
a_pj = delta\(Y_v'*M*x);
s_pj = Y_v*a_pj;
if min(s_pj)<=0
  warning('has negative or infinite conformal factors')
  
  az = 0:0.1:2*pi;
  el = 0:0.1:pi;
  [AZ,EL] = meshgrid(az,el);
  R = ones(size(AZ));
  [vx,vy,vz] = sph2cart(AZ(:),EL(:)-pi/2,R(:));
  vv = [vx, vy, vz];
  Y = sphericalHarmonicBase(vv,maxL);
  shsh = Y*a_pj;

  min(shsh)
end

%%
if strcmp('log',discriptor)
  modifier = @(x) exp(x);
elseif strcmp('inv',discriptor)
  modifier = @(x) log(abs(1./x))*(numel(find(x<0,1)) == 1)+ log((1./x))*isempty(find(x<0,1));
else
  modifier = @(x) x;
end

crange = [min(modifier(s_pj)) max(modifier(s_pj))];
figure(); colormap jet
set(gcf,'outerposition',[0, 0, 1920, 1080]);

subplot(2,2,1); hold all; axis equal; view([0 90]);
title('cMCF conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_T),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('eastoutside');
ylabel(ch,'s');

subplot(2,2,2); hold all; axis equal; view([0 90]);
title('cMCF mesh grid (hammered)');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3));
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('westoutside');
ylabel(ch,'s');

subplot(2,2,3); hold all; axis equal; view([0 90]);
title('SH projected conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_pj),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('eastoutside');
ylabel(ch,'s');

subplot(2,2,4); hold all; axis equal; view([0 90]);
title('SH least-squares fitted conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_ls),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('westoutside');
ylabel(ch,'s');
%%
% hgexport(gcf,num2str(maxL,[target 'SphericalHarmonicsL=%02d' discriptor '.png']),...
%   hgexport('factorystyle'), 'Format', 'png');
pause(.1);
apjlist = padcat(apjlist,a_pj);
alslist = padcat(alslist,a_ls);
end
% save(['SH/forward/' target 'SphericalHarmonics.' discriptor '.mat'],'apjlist','alslist');

%% add caption and make gif
% unix(['mogrify -font Liberation-Sans -fill white -undercolor ''#000000F0'' -pointsize 26 ' ...
%   '-gravity NorthEast -annotate +10+10 %t ' target '*phericalH*' discriptor '.png']);
% 
% unix(['convert -delay 100 -loop 0 ' ...
%   target '*phericalH*' discriptor '.png SH/forward/' target 'SphericalHarmonics.' discriptor '.gif']);
% 
% unix(['rm ' target '*phericalH*' discriptor '.png']);

%% test if spherical harmonics look like spherical harmonics 
%  (this debugging procedure failed me because of all the symmetries of SH)
% a = zeros((maxL+1)^2,1);
% a(5) = 2*sqrt(pi);
% sh = Y_v*a;
% 
% figure();
% trisurf(f_T,v(:,1),v(:,2),v(:,3),sh,'facecolor','interp','edgecolor','none'); axis equal;
% colorbar; colormap jet

%% get original spectrum
numeig = max(maxL+1)^2;
[M_T,L_T] = lapbel(v_T,f_T);
D_T = eigvf(L_T,M_T,numeig);

% [M,L] = lapbel(v,f_T);
% D_w = eigvf(L,diag(sparse(s_T)).*M,numeig);

%% get SH basis spectrum
ai = 1;
for maxL = Lrange
numeig = max(maxL+1)^2;
a_pj = apjlist(1:numeig,ai);
a_ls = alslist(1:numeig,ai);
ai = ai + 1;
D_s = [];
for l = 0:maxL
    D_s = [D_s, repmat(-l*(l+1),1,2*l+1)]; % spherical harmonic eigenvalues
end

L_pj = zeros(numeig);
L_ls = zeros(numeig);
for i = 1:numeig
  if exist(num2str(i,'../RSHI/RSHI%04d.mat'),'file')
    load(num2str(i,'../RSHI/RSHI%04d.mat'));
    computeNow = 0;
  else
    warning('ain''t nobody got time for this');
    computeNow = 1;
  end
  if computeNow
    for j = 1:numeig
      apjkcijk = 0;
      alskcijk = 0;
      for k = 1:numeig
        cijk = integrate3realSH(LM([i j k],:));
        apjkcijk = apjkcijk + a_pj(k)*cijk;
        alskcijk = alskcijk + a_ls(k)*cijk;
      end
      L_pj(i,j) = D_s(i)*apjkcijk;
      L_ls(i,j) = D_s(i)*alskcijk;
    end
  else
    L_pj(i,:) = D_s(i)*a_pj'*real(cijk(1:numeig,1:numeig));
    L_ls(i,:) = D_s(i)*a_ls'*real(cijk(1:numeig,1:numeig));
  end
end

D_pj = sort(eig(L_pj));
D_ls = sort(eig(L_ls));

%% compare spectrum
rei = numeig:-1:1;
close all;figure();hold all;
plot(rei,-D_pj,'--','linewidth',2)
plot(rei,-D_ls,':','linewidth',2)

plot(rei,-D_T(end-numeig+1:end))
if strfind(target,'spot')
  load spotspec.mat
elseif strfind(target,'bunny')
  load bunnyspec.mat
else
  warning('no discrete spectrum data available');
  D = [];
end
if size(D,1) > numel(rei)
  rei = [nan(size(size(D,1):-1:numeig+1)) rei];
else
  rei = rei(end-size(D,1)+1:end);
end
for i = 1:size(D,2)
  plot(rei,-D(:,i));
end
xlabel('# of Eigenvalues');
ylabel('Laplacian Eigenvalues');
title(num2str(maxL,'Conformal Factor Forward Problem in Spherical Harmonics L=%02d'));
legend('projected SH spectrum','least squares SH spectrum','discrete cotan spectrum',...
  'location','best');
% saveas(gcf,num2str(maxL,['SH/forward/' target 'SHspecL=%d.png']));
%%
% close all;
figure(); hold all; grid on;
plot(numeig:-1:1,abs((-D_pj+D_T)./D_T));
% plot(numeig:-1:1,abs((-D_ls+D_T)./D_T));
set(gca,'yscale','log');
xlabel('# of Eigenvalues');
ylabel('log Laplacian Eigenvalues');
saveas(gcf,num2str(maxL,['SH/forward/' target 'SHspecLogDiffL=%d.png']));

save(num2str(maxL,['SH/forward/' target 'SHspecL=%d.mat']),'a_pj','D_pj');
end