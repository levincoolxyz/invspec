clear;

% target = 'spot';
target = 'bunny';

% discriptor = 'inv';
discriptor = 'log';

if exist(['mcf/' target '.mat'],'file')
  load(['mcf/' target '.mat'])
else
  fid = fopen(['../meshes/' target '.obj'],'rt');
  [v_T,f_T] = readwfobj(fid);
  [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');
  save(['mcf/' target '.mat'],'v_T','f_T','s_T','v');
end

[vham,fham] = xyz2hammer(v);
%%
apjlist = [];
alslist = [];
for L = [1:3 5 8 10 15 20 25]

Y_v = sphericalHarmonicBase(v,L);

if strcmp('log',discriptor)
  x = log(1./s_T);
else
  x = s_T;
end

a_ls = (Y_v'*Y_v)\(Y_v'*x);
s_ls = Y_v*a_ls;

M = lapbel(v,f_T);
idtest = Y_v'*M*Y_v;
a_pj = Y_v'*M*x;
s_pj = Y_v*a_pj;

if strcmp('log',discriptor)
  s_ls = 1./exp(s_ls);
  s_pj = 1./exp(s_pj);
end

%%
% modifier = @(x) log(1./x);
modifier = @(x) x;
close all; figure(); colormap jet
set(gcf,'outerposition',[0, 0, 1920, 1080]);
subplot(2,2,1); hold all; axis equal; view([0 90]);
title('cMCF conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_T),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');
subplot(2,2,2); hold all; axis equal; view([0 90]);
title('cMCF mesh grid (hammered)');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3));
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');
subplot(2,2,3); hold all; axis equal; view([0 90]);
title('SH projected conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_pj),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');
subplot(2,2,4); hold all; axis equal; view([0 90]);
title('SH least-squares fitted conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_ls),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

%%
hgexport(gcf,num2str(L,[target 'SphericalHarmonicsL=%02d' discriptor '.png']),...
  hgexport('factorystyle'), 'Format', 'png');
apjlist = padcat(apjlist,a_pj);
alslist = padcat(alslist,a_ls);
end
save([target 'SphericalHarmonics.' discriptor '.mat'],'apjlist','alslist','Y_v');

%% test if spherical harmonics look like spherical harmonics
% a = zeros(size(a_pj));
% a(18) = 1;
% 
% sh = Y_v*a;
% figure();
% trisurf(f_T,v(:,1),v(:,2),v(:,3),sh,'facecolor','interp','edgecolor','none'); axis equal;

%% add caption and make gif
unix(['mogrify -font Liberation-Sans -fill white -undercolor ''#000000F0'' -pointsize 26 ' ...
  '-gravity NorthEast -annotate +10+10 %t ' target '*phericalH*' discriptor '.png']);

unix(['convert -delay 100 -loop 0 ' ...
  target '*phericalH*' discriptor '.png ' target 'SphericalHarmonics.' discriptor '.gif']);

unix(['rm ' target '*phericalH*' discriptor '.png']);