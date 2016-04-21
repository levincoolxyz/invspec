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
  colormap('jet');
  saveas(figh,[token '.png']);
%   [token,remain] = strtok(remain);
  close all;
end

%% two gifs
unix('convert -delay 100 -loop 0 i2_200_t2_abs\(Y32\(v\)\)_e*.png i2_300_t2_abs\(Y32\(v\)\)_e0.6p2.png i2_200-300_t2_abs\(Y32\(v\)\)_e0.6p0.4-2.gif');
unix('convert -delay 100 -loop 0 i3_300_t2_abs\(Y33\(v\)\)_e*.png i3_300_t2_abs\(Y33\(v\)\)_e0.1-1p0.512.gif');
