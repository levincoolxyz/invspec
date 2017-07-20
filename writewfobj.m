function fid = writewfobj(meshname,v,f)
fid = fopen(['../meshes/' meshname '.obj'],'w');
fprintf(fid,'v %f %f %f\n',v');
fid = fopen(['../meshes/' meshname '.obj'],'a');
fprintf(fid,'f %d// %d// %d//\n',f');
end