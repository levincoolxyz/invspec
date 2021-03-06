function [figh] = visualizeSH(v,v_T,v_end,vmcf,f,f_T,a_end,s_end,s_T,...
  J_hist,Jc_hist,D_0,D_T,D_endp,D_end,pcable)
% function [figh] = visualizeSH(v,v_T,v_end,f,f_T,s_end,s_T,...
%   J_hist,Jc_hist,D_0,D_T,D_endp,D_end,pcable)
% 
% pcable - true (1) if want to use principal component analysis to align meshes

if nargin<14 || isempty(pcable), pcable = 1; end
numv = size(v,1);
%% principal component analysis for alignment
if pcable
  dcm = pa(v,f); % princomp(v); %pca(v);
  if ~isempty(v_T)
    dcm_T = pa(v_T,f_T); % princomp(v_T); %pca(v_T);
  end
  dcm_end = pa(v_end,f); % princomp(v_end); %pca(v_end);

  % ensure orientation preserving
  dcm = dcm*det(dcm);
  if ~isempty(v_T), dcm_T = dcm_T*det(dcm_T); end
  dcm_end = dcm_end*det(dcm_end);

  % rotate vertices
  v = dcmrot(v,dcm);
  if ~isempty(v_T), v_T = dcmrot(v_T,dcm_T); end
  v_end = dcmrot(v_end,dcm_end);
end
%% show SH coefficient decay
if ~isempty(a_end)
  a = nonzeros(a_end);
  L = 0:sqrt(numel(a))-1;
  amean = zeros(size(L));
  aup = zeros(size(L));
  adown = zeros(size(L));
  for l = L
    al = (abs(a((1:2*l+1)+l^2)));
    amean(l+1) = mean(al);
    aup(l+1) = max(al)-mean(al);
    adown(l+1) = mean(al)-min(al);
  end
  figure(); hold all; grid on;
  errorbar(L,log10(amean),.434*(adown./amean),.434*(aup./amean),'r:x')
  plot(L,log10(1./L),'k--');
  plot(L,log10(1./L.^2),'k--');
  legend('mean','L^{-1}','L^{-2}');
  set(gca,'xlim',[0 max(L)]);
  ylabel('Coefficient Magnitude $log(|a|)$','interpreter','latex');
  xlabel('SH Degree $L$','interpreter','latex');
  title('Spherical Harmonic Coefficients Decay vs. its Degree');
end
%% initialize
figh = figure();
set(gcf,'outerposition',[0, 0, 1920, 1080]);
vlim = [];
colormap('jet');
%% compare mesh
if isempty(v_T)
  subplot(2,4,1);
  title('Target spectrum');
  set(gca, 'visible', 'off')
else
  subplot(2,4,1); hold all; grid on; axis equal
%   view(3); 
  view([-68 68]);set(gca,'fontsize',16)
  trimesh(f_T,v_T(:,1),v_T(:,2),v_T(:,3));
  vlim = max(max(abs(v_T)));
  vlim = [-vlim vlim];
  set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
  xlabel('x'); ylabel('y'); zlabel('z');
  title('Target mesh');
end

subplot(2,4,4); hold all; grid on; axis equal
%   view(3); 
  view([-68 68]);set(gca,'fontsize',16)
if ~isempty(v_end)
  trimesh(f,v_end(:,1),v_end(:,2),v_end(:,3));
  vlim = max([vlim max(max(abs(v_end)))]);
  vlim = [-vlim vlim];
  set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
  xlabel('x'); ylabel('y'); zlabel('z');
else
  text(0,0,0,'N/A'); grid off; axis off;
end
title('Resultant mesh');

%% compare conformal factors
vlim = max(max(abs(v)));
vlim = [-vlim vlim];

if isempty(vmcf)
  [s_T,vmcf] = meancurvflow(v_T,f_T,1e5,'c');
end

% crange = [min([s_T;s_end]) max([s_T;s_end])];
crange = [min([s_T;]) max([s_T;])];
subplot(2,4,2); hold all; grid on; axis equal
title('cMCF conformal factors');
set(gca,'fontsize',16)

% 3d spherical plot
% trisurf(f_T,v_T_mcf(:,1),v_T_mcf(:,2),v_T_mcf(:,3),s_T,...
%   'facecolor','interp','edgecolor','none');
% set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
% xlabel('x'); ylabel('y'); zlabel('z'); view(3); grid on; 

% hammer projection plot
[vham1,f_T_ham] = xyz2hammer(vmcf);
trisurf(f_T_ham,vham1(:,1),vham1(:,2),vham1(:,3),s_T,...
  'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on'); % to see title etc.

caxis(crange);
ch = colorbar('southoutside');
ylabel(ch,'1/s');

subplot(2,4,3); hold all; grid on; axis equal
title('Spectrally optimized factors');
set(gca,'fontsize',16)

% 3d spherical plot
% trisurf(f,v(:,1),v(:,2),v(:,3),s_end,...
%   'facecolor','interp','edgecolor','none');
% set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
% xlabel('x'); ylabel('y'); zlabel('z'); view(3);

% hammer projection plot
[vham2,f_ham] = xyz2hammer(v);
trisurf(f_ham,vham2(:,1),vham2(:,2),vham2(:,3),s_end,...
  'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on'); % to see title etc.

caxis(crange);
ch = colorbar('southoutside');
ylabel(ch,'1/s');
%% compare spectra
numeig = numel(D_T);
% [M_T,L_T] = lapbel(v_T,f_T);
% D_T = eigvf(L_T,M_T,numeig);
% [M,L] = lapbel(v,f);
% D0 = eigvf(L,M,size(v,1));
[Mmcf,Lmcf] = lapbel(vmcf,f_T);
D_MCF = eigvf(Lmcf,sparse(diag(1./s_T))*Mmcf,numeig);

v0 = (D_MCF - D_T)./D_T;
v1 = (D_endp - D_T)./D_T;
if ~isempty(D_end)
  v2 = (D_end - D_T)./D_T;
else
  v2 = [];
end
% v3 = (D0 - D_T)./D_T;

% reverse eigenvalue indexing (change neg-def L convention to pos-def L)
rei = size(v1,1):-1:2;

subplot(2,4,5:8); hold all; grid on;
% linear relative error plot
% plot(rei, v0(1:end-1),'b--','linewidth',2);
% plot(rei, v1(1:end-1),'k-','linewidth',2);
% plot(rei, v2(1:end-1),'ro');
% plot(rei, v3(1:end-1),'g--','linewidth',2);
% c1 = .18;
% c2 = 0;
% legendpos = 'northeast';

% log |relative error| plot
plot(rei, abs(v0(1:end-1)),'b--','linewidth',2);
plot(rei, abs(v1(1:end-1)),'k-','linewidth',2);
if ~isempty(v2), plot(rei, abs(v2(1:end-1)),'ro'); end
% plot(rei, abs(v3(1:end-1)),'g--','linewidth',2);
set(gca,'yscale','log');
c1 = 0;
c2 = 0.18;
% legendpos = 'southeast';
legendpos = 'best';

legend('(\lambda_{cMCF} - \lambda_{target})/\lambda_{target}',...
  '(\lambda_{MIEP2} - \lambda_{target})/\lambda_{target}',...
  '(\lambda_{final embed} - \lambda_{target})/\lambda_{target}',...
  'location',legendpos);
xlabel('# of eigenvalues (in ascending magnitudes)');
title('Deviation from target Laplacian eigenvalues');

set(gca,'xlim',[2 numel(D_T)])
xm = get(gca,'xlim');
ym = get(gca,'ylim');
ym = [ym(1) 2]; set(gca,'ylim',ym);
set(gca,'fontsize',16)
set(gca,'xscale','log');
% text(floor(max(xm)/4),((ym(2)/ym(1))^c2*ym(2) + c1*diff(ym)),...
text(floor(max(xm)^.4),((ym(2)/ym(1))^c2*ym(2) + c1*diff(ym)),...
  num2str([J_hist(end) Jc_hist(end) numeig],...
  ['J_{MIEP2} = %g     ',...
  'J_{embedding} = %g   numeig=%d']));
end