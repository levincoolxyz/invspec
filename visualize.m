function [figh] = visualize(v,v_T,v_end,f,f_T,s_end,s_T,...
  J_hist,Jc_hist,D_0,D_T,D_endp,D_end,pcable)
% function [figh] = visualize(v,v_T,v_end,f,f_T,s_end,s_T,...
%   J_hist,Jc_hist,D_0,D_T,D_endp,D_end,pcable)
% 
% pcable - true (1) if want to use principal component analysis to align meshes

if nargin<14 || isempty(pcable), pcable = 1; end
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
%% initialize
figh = figure();
set(gcf,'outerposition',[0, 0, 1920, 1080]);
vlim = [];
colormap('jet');
crange = [min([s_T;s_end]) max([s_T;s_end])];
%% compare mesh
if isempty(v_T)
  subplot(2,4,1);
  title('Target spectrum');
  set(gca, 'visible', 'off')
else
  v_T = v_T - mean(v_T);
  subplot(2,4,1); hold all; grid on; axis equal
  view(3); 
  view([54 54]);
%   view([-68 68]);
  set(gca,'fontsize',16)
%   trimesh(f_T,v_T(:,1),v_T(:,2),v_T(:,3));
  trisurf(f_T,v_T(:,1),v_T(:,2),v_T(:,3),s_T);
  vlim = max(max(abs(v_T)));
  vlim = [-vlim vlim];
  set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
  xlabel('x'); ylabel('y'); zlabel('z');
  title('Target mesh');
end

v_end = v_end - mean(v_end);
subplot(2,4,4); hold all; grid on; axis equal
  view(3); 
  view([54 54]);
%   view([-68 68]);
  set(gca,'fontsize',16)
% trimesh(f,v_end(:,1),v_end(:,2),v_end(:,3));
trisurf(f,v_end(:,1),v_end(:,2),v_end(:,3),s_end);
vlim = max([vlim max(max(abs(v_end)))]);
vlim = [-vlim vlim];
set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
xlabel('x'); ylabel('y'); zlabel('z');
title('Resultant mesh');

%% compare conformal factors
% vlim = max(max(abs(v)));
% vlim = [-vlim vlim];

if numel(s_T) ~= size(v_T,1)
  [s_T,v_T_mcf] = meancurvflow(v_T,f_T,1e5,'c');
else
  v_T_mcf = v;
end

subplot(2,4,2); hold all; grid on; axis equal
title('cMCF conformal factors');
set(gca,'fontsize',16)

% 3d spherical plot
% trisurf(f_T,v_T_mcf(:,1),v_T_mcf(:,2),v_T_mcf(:,3),s_T,...
%   'facecolor','interp','edgecolor','none');
% set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
% xlabel('x'); ylabel('y'); zlabel('z'); view(3); grid on; 

% hammer projection plot
[vham1,f_T_ham] = xyz2hammer(v_T_mcf);
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
[M,L] = lapbel(v,f);
D0 = eigvf(L,M,size(v,1));
if size(s_T,1) ~= size(M,1)
  s_T = ones(size(M,1),1);
end
D_MCF = eigvf(L,diag(1./s_T)*M,size(v,1));
[M_T,L_T] = lapbel(v_T,f_T);
D_T = eigvf(L_T,M_T,size(v,1));
[M_end,L_end] = lapbel(v_end,f);
D_end = eigvf(L_end,M_end,size(v,1));
if ~isempty(s_end)
  D_endp = eigvf(L,diag(1./s_end)*M,size(v,1));
else
  D_endp = [nan(size(v,1)-numel(D_endp),1); D_endp];
end

subplot(2,4,5:8); hold all; grid on;
if size(M,1) >= numel(D_T)
  D_T = [nan(size(M,1)-numel(D_T),1);D_T];
elseif size(M,1) <= numel(D_T)
  D_T = D_T(end-size(M,1)+1:end);
end
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

% log |relative error| plot,...
plot(rei, abs(v3(1:end-1)),'g--','linewidth',2,'displayname',...
  '|1-\lambda_{init}/\lambda_{T}|');
plot(rei, abs(v0(1:end-1)),'b--','linewidth',2,'displayname',...
  '|1-\lambda_{cMCF}/\lambda_{T}|');
plot(rei, abs(v1(1:end-1)),'k-','linewidth',2,'displayname',...
  '|1-\lambda_{spec}/\lambda_{T}|');
plot(rei, abs(v2(1:end-1)),'ro--','displayname',...
  '|1-\lambda_{embed}/\lambda_{T}|');
set(gca,'yscale','log');
c1 = 0;
c2 = 0.18;
legendpos = 'best';

legend('location',legendpos);
xlabel('# of eigenvalues (in ascending magnitudes)');
title('Relative deviation from target spectrum \lambda_{T}');

set(gca,'xlim',[2 numel(D_T)])
xm = get(gca,'xlim');
ym = get(gca,'ylim');
ym = [ym(1) 2]; set(gca,'ylim',ym);
set(gca,'fontsize',16)
set(gca,'xscale','log');
eigencost(log(s_end),M,L,D_T,numeig,0.01)
% text(floor(max(xm)/45),((ym(2)/ym(1))^c2*ym(2) + c1*diff(ym)),...
text(floor(max(xm)^.45),((ym(2)/ym(1))^c2*ym(2) + c1*diff(ym)),...
  num2str([.5*sum(v1(end-numeig:end-1).^2) Jc_hist(end) numeig],...
  ['J_{spec} = %g     ',...
  'J_{embed} = %g   numeig=%d']));
end