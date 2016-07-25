clear;
%% load target and mcf data
% target = 'spot';
% target = 'bunny';
target = 'spotS10';

discriptor = '';
% discriptor = 'inv';
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
for L = [1:3 5 8 10 15 20 25]
%% get spherical harmonics on vertices
[Y_v,LM] = sphericalHarmonicBase(v,L);
numSH = size(LM,1);

%% calculate conf. factors in SH basis
x = s_T;

a_ls = (Y_v'*Y_v)\(Y_v'*x);
s_ls = Y_v*a_ls;

M = lapbel(v,f_T);
idtest = Y_v'*M*Y_v;
a_pj = Y_v'*M*x;
s_pj = Y_v*a_pj;

%%
if strcmp('log',discriptor)
  modifier = @(x) exp(x);
elseif strcmp('inv',discriptor)
  modifier = @(x) 1./x;
else
  modifier = @(x) x;
end

crange = [min(modifier(s_pj)) max(modifier(s_pj))];
close all; figure(); colormap jet
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
hgexport(gcf,num2str(L,[target 'SphericalHarmonicsL=%02d' discriptor '.png']),...
  hgexport('factorystyle'), 'Format', 'png');
apjlist = padcat(apjlist,a_pj);
alslist = padcat(alslist,a_ls);
end
save(['SH/' target 'SphericalHarmonics.' discriptor '.mat'],'apjlist','alslist');

%% add caption and make gif
unix(['mogrify -font Liberation-Sans -fill white -undercolor ''#000000F0'' -pointsize 26 ' ...
  '-gravity NorthEast -annotate +10+10 %t ' target '*phericalH*' discriptor '.png']);

unix(['convert -delay 100 -loop 0 ' ...
  target '*phericalH*' discriptor '.png SH/' target 'SphericalHarmonics.' discriptor '.gif']);

unix(['rm ' target '*phericalH*' discriptor '.png']);

%% test if spherical harmonics look like spherical harmonics
% a = zeros(size(a_pj));
% a(18) = 2*sqrt(pi);

% sh = Y_v*a;
% figure();
% trisurf(f_T,v(:,1),v(:,2),v(:,3),sh,'facecolor','interp','edgecolor','none'); axis equal;

%% get original spectrum
% numeig = max(L+1)^2;
% [M_T,L_T] = lapbel(v_T,f_T);
% D_T = eigvf(L_T,M_T,numeig);
% 
% % [M,L] = lapbel(v,f_T);
% % D_w = eigvf(L,diag(sparse(s_T)).*M,numeig);
% 
% %% get SH basis spectrum
% % a = zeros(size(a_pj));
% % a(1) = 2*sqrt(pi);
% % aa = a;
% % aa = [2*sqrt(pi), -1.0066359, 0.57527065, 0.31905909, -0.17029902, ...
% % -0.00090668897, 0.0058173331, 0.0031209552, 0.0048153206, ...
% % 0.0038694581, 0.00016732501, 0.00030969496, -0.0002770344, ...
% % 0.000070761836, -0.00020614028, -0.000039104542];
% aa = a_pj;
% % aa = a_ls;
% 
% D_s = [];
% for l = 0:L
%     D_s = [D_s, repmat(-l*(l+1),1,2*l+1)]; % spherical harmonic eigenvalues
% end
% 
% L_sh = zeros(numeig);
% for i = 1:numeig
%   for j = 1:numeig
%     akcijk = 0;
%     for k = 1:numeig
%       cijk = integrate3realSH(LM([i j k],:));
%       akcijk = akcijk + aa(k)*cijk;
%     end
%     L_sh(i,j) = D_s(i)*akcijk;
%   end
% end
% 
% D_sh = eig(L_sh);
% if norm(imag(D_sh)) > 1e-10
%   error('no no no no no no ...');
% end
% D_sh = sort(real(D_sh));
% 
% %% compare spectrum
% rei = numeig:-1:1;
% figure();hold all;
% plot(rei,-D_T)
% % plot(rei,-D_w,'kx')
% plot(rei,-D_sh,'--')
% 
% %% more comparisons
% % % load cowspec.mat
% % load bunspec.mat
% % rei = [nan(size(size(D,1):-1:numeig+1)) rei];
% % for i = 1:size(D,2)
% %   plot(rei,-D(:,i));
% % end