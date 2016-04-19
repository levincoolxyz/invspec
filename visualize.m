function [figh] = visualize(v,f,v_end,v_T,f_T,...
  D_0,D_T,D_endp,D_end)
% compare mesh
Mesh0 = TriRep(f,v); %triangulation(f,v);
Mesh_T = TriRep(f_T,v_T); %triangulation(f_T,v_T);
Mesh_end = TriRep(f,v_end); %triangulation(f,v_end);

figh = figure(); set(gcf,'outerposition',[0, 0, 1024, 768]);
subplot(2,3,1); hold all; view(3); grid on; axis equal
trimesh(Mesh0);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('original mesh');

subplot(2,3,2); hold all; view(3); grid on; axis equal
trimesh(Mesh_T);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('target mesh');

subplot(2,3,3); hold all; view(3); grid on; axis equal
trimesh(Mesh_end);
set(gca,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
title('resultant mesh');
view([-140 20]);

% compare spectra
subplot(2,3,4:6); hold all; grid on;
plot((D_0 - D_T),'bx:');
% plot((D_Tp - D_T),'gs');
plot((D_endp - D_T),'k-','linewidth',2);
plot((D_end - D_T),'ro');
% if exist('D_Tp','var')
%   legend('\lambda_{initial} - \lambda_{target}',...
%     '\lambda_{target embed} - \lambda_{target}',...
%     '\lambda_{MIEP2} - \lambda_{target}',...
%     '\lambda_{final embed} - \lambda_{target}',...
%     'location','best');
% else
  legend('\lambda_{initial} - \lambda_{target}',...
    '\lambda_{MIEP2} - \lambda_{target}',...
    '\lambda_{final embed} - \lambda_{target}',...
    'location','best');
% end
% set(gca,'yscale','log');
xlabel('# of eigenvalues (#1 is of the highest frequency)'); 
ylabel('\propto eigenvalue magnitude');
title('Deviation from target Laplacian eigenvalues in magnitude');

ym = get(gca,'ylim');
xm = get(gca,'xlim');
text(floor(max(xm)/4.5),max(ym(2) + .18*diff(ym)),...
  num2str([J_hist(end) Jc_hist(end)],...
  ['Convergence Energies: J_{MIEP2} = %g     ',...
  'J_{embedding} = %g']));

end