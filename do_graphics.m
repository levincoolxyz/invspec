close all;
folder = 'bunny/recursive602';
[~,datalist] = unix(['basename -s .mat -a ' folder '/i*.mat']);
% [token,remain] = strtok(datalist);
datalist = textscan(datalist,'%s');
for i = 1:size(datalist{1},1)
% while size(token,1)>=1
  s_end = []; s_T = []; % backward compatiability
  token = char(datalist{1}(i));
  load([folder '/' token '.mat']);
  figh = visualize(v,v_T,v_end,f,f_T,s_end,s_T,...
    J_hist,Jc_hist,D_0,D_T,D_endp,D_end,0);
  hgexport(figh,[folder '/' token '.png'],...
      hgexport('factorystyle'), 'Format', 'png');
%   [token,remain] = strtok(remain);
 close all;
end

%% two gifs
if strcmp(folder,'Y32')
unix('mv Y32/i2_300_t2_abs\(Y32\(v\)\)_e0.1p1.png Y32/i2_300_t2_abs\(Y32\(v\)\)_e0.1p1.0.png');
unix('mv Y32/i2_300_t2_abs\(Y32\(v\)\)_e0.1p2.png Y32/i2_300_t2_abs\(Y32\(v\)\)_e0.1p2.0.png');
unix('convert -delay 100 -loop 0 Y32/i2_300_t2_abs\(Y32\(v\)\)_e*.png i2_300_t2_abs\(Y32\(v\)\)_e0.1p0.5-2.gif');
unix('mv Y32/i3_300_t2_abs\(Y32\(v\)\)_e0.1p1.png Y32/i3_300_t2_abs\(Y32\(v\)\)_e0.1p1.0.png');
unix('mv Y32/i3_300_t2_abs\(Y32\(v\)\)_e0.1p2.png Y32/i3_300_t2_abs\(Y32\(v\)\)_e0.1p2.0.png');
unix('convert -delay 100 -loop 0 Y32/i3_300_t2_abs\(Y32\(v\)\)_e*.png i3_300_t2_abs\(Y32\(v\)\)_e0.1p0.5-2.gif');
elseif strcmp(folder,'Y33')
unix('mv Y33/i3_300_t2_abs\(Y33\(v\)\)_e1p0.5.png Y33/i3_300_t2_abs\(Y33\(v\)\)_e1.0p0.5.png');
unix('convert -delay 100 -loop 0 Y33/i3_300_t2_abs\(Y33\(v\)\)_e*.png i3_300_t2_abs\(Y33\(v\)\)_e0.1-1p0.5.gif');
unix('mv Y33/i2_300_t2_abs\(Y33\(v\)\)_e1p0.5.png Y33/i2_300_t2_abs\(Y33\(v\)\)_e1.0p0.5.png');
unix('convert -delay 100 -loop 0 Y33/i2_300_t2_abs\(Y33\(v\)\)_e*.png i2_300_t2_abs\(Y33\(v\)\)_e0.1-1p0.5.gif');
end