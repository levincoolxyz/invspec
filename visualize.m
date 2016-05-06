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
if ~isempty(v_T)
  dcm_T = dcm_T*det(dcm_T);
end
dcm_end = dcm_end*det(dcm_end);
% rotate vertices
v = dcmrot(v,dcm);
if ~isempty(v_T)
  v_T = dcmrot(v_T,dcm_T);
end
v_end = dcmrot(v_end,dcm_end);

figh = figure(); set(gcf,'outerposition',[0, 0, 1024, 768]);
%% compare mesh
Mesh0 = TriRep(f,v); %triangulation(f,v);
subplot(2,3,1); hold all; view(3); grid on; axis equal
trimesh(Mesh0);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('original mesh');

if isempty(v_T)
  subplot(2,3,2);
  text(0,0.5,'target given as spectrum');
  set(gca, 'visible', 'off')
else
  Mesh_T = TriRep(f_T,v_T); %triangulation(f_T,v_T);
  subplot(2,3,2); hold all; view(3); grid on; axis equal
  trimesh(Mesh_T);
  set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
  xlabel('x'); ylabel('y'); zlabel('z');
  title('target mesh');
end

Mesh_end = TriRep(f,v_end); %triangulation(f,v_end);
subplot(2,3,3); hold all; view(3); grid on; axis equal
trimesh(Mesh_end);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('resultant mesh');

%% compare spectra

subplot(2,3,4:6); hold all; grid on;
v1 = (D_endp - D_T)./D_T;
v2 = (D_end - D_T)./D_T;
plot(v1(1:end-1),'k-','linewidth',2);
plot(v2(1:end-1),'ro');
legend('(\lambda_{MIEP2} - \lambda_{target})/\lambda_{target}',...
  '(\lambda_{final embed} - \lambda_{target})/\lambda_{target}',...
  'location','northwest');
% plot((D_endp - D_T),'k+-','linewidth',2);
% plot((D_end - D_T),'ro:');
% legend('\lambda_{MIEP2} - \lambda_{target}',...
%   '\lambda_{final embed} - \lambda_{target}',...
%   'location','southeast');
xlabel('# of eigenvalues (#1 is of the highest frequency)');
% ylabel('\propto eigenvalue magnitude');
title('Deviation from target Laplacian eigenvalues');

fig = gcf;
if isstruct(fig)
  ax = fig.CurrentAxes;
else
  ax = gca;
end
pause(.1);
set(ax,'xlim',[1 numel(D_0)])
xm = get(ax,'xlim');
ym = get(ax,'ylim');
text(floor(max(xm)/4),max(ym(2) + .18*diff(ym)),...
  num2str([numel(J_hist) J_hist(end) numel(Jc_hist) Jc_hist(end)],...
  ['iter#%d: J_{MIEP2} = %g     ',...
  'iter#%d: J_{embedding} = %g']));

colormap('jet');

%% correlate conformal factors

end
