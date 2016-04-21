close all;
[~,datalist] = unix('basename -s .mat -a i*.mat');
% [token,remain] = strtok(datalist);
datalist = textscan(datalist,'%s');
for i = 1:size(datalist{1},1)
% while size(token,1)>=1
  token = char(datalist{1}(i));
  load([token '.mat']);
  figh = visualize(v,f,v_end,v_T,f_T,...
  J_hist,Jc_hist,D_0,D_T,D_endp,D_end);
  saveas(figh,[token '.png']);
%   [token,remain] = strtok(remain);
  close all;
end
