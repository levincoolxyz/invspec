close all;
load spot/recursive487/i4_mcf_t3_spot487_e0.95p0r0.1.mat

[vham,fham] = xyz2hammer(v);
alist = [];

for L = 3%[1:3 5 8 10 15 20 25 29]

[a,s_T_sh] = sphericalHarmonicDecomposition(v,s_T,L);
figure(); colormap jet

subplot(1,2,1); hold all; axis equal; view([0 90]);
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),s_T,'facecolor','interp','edgecolor','none');
set(gca,'visible','off');
subplot(1,2,2); hold all; axis equal; view([0 90]);
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),s_T_sh,'facecolor','interp','edgecolor','none');
set(gca,'visible','off');

% % hgexport(gcf,num2str(L,'bunSphericalHarmonicsL=%d.png'),...
% hgexport(gcf,num2str(L,'spotSphericalHarmonicsL=%d.png'),...
%   hgexport('factorystyle'), 'Format', 'png');
alist = padcat(alist,a);
end