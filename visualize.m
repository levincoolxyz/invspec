function [figh] = visualize(v,v_T,v_end,f,f_T,s_end,s_T,...
  J_hist,Jc_hist,D_0,D_T,D_endp,D_end)
% function [figh] = visualize(v,v_T,v_end,f,f_T,s_end,s_T,...
%   J_hist,Jc_hist,D_0,D_T,D_endp,D_end)

%% principal component analysis for alignment
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

%% initialize
figh = figure();
set(gcf,'outerposition',[0, 0, 1920, 1080]);
vlim = [];
colormap('jet');
%% compare mesh
if isempty(v_T)
  subplot(2,4,1);
  title('target spectrum');
  set(gca, 'visible', 'off')
else
  subplot(2,4,1); hold all; view(3); grid on; axis equal
  trimesh(f_T,v_T(:,1),v_T(:,2),v_T(:,3));
  vlim = max(max(abs(v_T)));
  vlim = [-vlim vlim];
  set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
  xlabel('x'); ylabel('y'); zlabel('z');
  title('target mesh');
end

subplot(2,4,4); hold all; view(3); grid on; axis equal
trimesh(f,v_end(:,1),v_end(:,2),v_end(:,3));
vlim = max([vlim max(max(abs(v_end)))]);
vlim = [-vlim vlim];
set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
xlabel('x'); ylabel('y'); zlabel('z');
title('resultant mesh');

%% compare conformal factors
crange = [min([s_T;s_end]) max([s_T;s_end])];
subplot(2,4,2); hold all; view(3); grid on; axis equal
trisurf(f,v(:,1),v(:,2),v(:,3),s_T,...
  'facecolor','interp','edgecolor','none');
vlim = max(max(abs(v)));
vlim = [-vlim vlim];
set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
xlabel('x'); ylabel('y'); zlabel('z');
title('spherical/cMCF mesh');
caxis(crange);
colorbar('southoutside')

subplot(2,4,3); hold all; view(3); grid on; axis equal
trisurf(f,v(:,1),v(:,2),v(:,3),s_end,...
  'facecolor','interp','edgecolor','none');
vlim = max(max(abs(v)));
vlim = [-vlim vlim];
set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
xlabel('x'); ylabel('y'); zlabel('z');
title('spectrally opt. mesh');
caxis(crange);
colorbar('southoutside')
%% compare spectra
numeig = numel(D_T);
[M,L] = lapbel(v,f);
D0 = eigvf(L,M,size(v,1));
D_MCF = eigvf(L,diag(1./s_T)*M,size(v,1));
[M_T,L_T] = lapbel(v_T,f_T);
D_T = eigvf(L_T,M_T,size(v,1));
[M_end,L_end] = lapbel(v_end,f);
D_end = eigvf(L_end,M_end,size(v,1));
if ~isempty(s_end)
  [M,L] = lapbel(v,f);
  D_endp = eigvf(L,diag(1./s_end)*M,size(v,1));
else
  D_endp = [nan(size(v,1)-numel(D_endp),1); D_endp];
end

subplot(2,4,5:8); hold all; grid on;
v0 = (D_MCF - D_T)./D_T;
v1 = (D_endp - D_T)./D_T;
v2 = (D_end - D_T)./D_T;
v3 = (D0 - D_T)./D_T;

% reverse eigenvalue indexing (change neg-def L convention to pos-def L)
rei = size(v0,1):-1:2;

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
plot(rei, abs(v2(1:end-1)),'ro');
plot(rei, abs(v3(1:end-1)),'g--','linewidth',2);
set(gca,'yscale','log');
c1 = 0;
c2 = 0.18;
legendpos = 'southeast';

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
set(gca,'xscale','log');
% text(floor(max(xm)/4),((ym(2)/ym(1))^c2*ym(2) + c1*diff(ym)),...
text(floor(max(xm)^.4),((ym(2)/ym(1))^c2*ym(2) + c1*diff(ym)),...
  num2str([numel(J_hist) J_hist(end) numel(Jc_hist) Jc_hist(end), ...
  numeig],...
  ['iter#%d: J_{MIEP2} = %g     ',...
  'iter#%d: J_{embedding} = %g   numeig=%d']));
end
