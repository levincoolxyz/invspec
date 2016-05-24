close all;
[~,datalist] = unix('basename -s .mat -a i*bunny2k*.mat');
% [token,remain] = strtok(datalist);
datalist = textscan(datalist,'%s');
for i = 1:size(datalist{1},1)
% while size(token,1)>=1
  s_end = []; s_T = []; % backward compatiability
  token = char(datalist{1}(i));
  load([token '.mat']);
  figh = visualize(v,v_T,v_end,f,f_T,s_end,s_T,...
    J_hist,Jc_hist,D_0,D_T,D_endp,D_end);
  pause(.1);
  hgexport(figh,[token '.png'],...
      hgexport('factorystyle'), 'Format', 'png');
%   [token,remain] = strtok(remain);
%   close all;
end

%% two gifs
% unix('mv i2_300_t2_abs\(Y32\(v\)\)_e0.1p1.png i2_300_t2_abs\(Y32\(v\)\)_e0.1p1.0.png')
% unix('mv i2_300_t2_abs\(Y32\(v\)\)_e0.1p2.png i2_300_t2_abs\(Y32\(v\)\)_e0.1p2.0.png')
% unix('convert -delay 100 -loop 0 i2_300_t2_abs\(Y32\(v\)\)_e*.png i2_300_t2_abs\(Y32\(v\)\)_e0.1p0.5-2.gif');
% unix('mv i3_300_t2_abs\(Y32\(v\)\)_e0.1p1.png i3_300_t2_abs\(Y32\(v\)\)_e0.1p1.0.png')
% unix('mv i3_300_t2_abs\(Y32\(v\)\)_e0.1p2.png i3_300_t2_abs\(Y32\(v\)\)_e0.1p2.0.png')
% unix('convert -delay 100 -loop 0 i3_300_t2_abs\(Y32\(v\)\)_e*.png i3_300_t2_abs\(Y32\(v\)\)_e0.1p0.5-2.gif');
% unix('mv i2_300_t2_abs\(Y33\(v\)\)_e1p0.5.png i2_300_t2_abs\(Y33\(v\)\)_e1.0p0.5.png')
% unix('convert -delay 100 -loop 0 i2_300_t2_abs\(Y33\(v\)\)_e*.png i2_300_t2_abs\(Y33\(v\)\)_e0.1-1p0.5.gif');
