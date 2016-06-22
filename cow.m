close all; clear
%% generate sequence of increasingly spherical cow
[v_T,f_T] = readwfobj(fopen('../meshes/spot487.obj','rt'));
[~,L_T0] = lapbel(v_T,f_T);
figure(); view(3); grid on; axis equal
for i = 1:40
  [~,v_T,~] = meancurvflow(v_T,f_T,1e-3,'c',L_T0,[],5);
  trimesh(f_T,v_T(:,1),v_T(:,2),v_T(:,3));
  save(num2str(i,'../meshes/cow%02d.mat'),'f_T','v_T');
  pause(.1);
end