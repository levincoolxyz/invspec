function [v,f] = readwfobj(fid)
% function [v,f] = readwfobj(fid)
% read wavefront object (*.obj) file
% 
% INPUT
% fid   - MATLAB file id pointing to the object file
% OUPUTS
% v     - (#vertices x 3) vertex matrix
% f     - (#faces x 3) face connectivity matrix
% 

temp = textscan(fid,'%s','Delimiter','');
fclose(fid);
temp = regexprep(temp{1,1},'\/[0-9]*\/?[0-9]*','');
v = regexp(temp,'^v .*','match');
f = regexp(temp,'^f .*','match');
clear temp;

v = strcat(v{~cellfun(@isempty,v)});
f = strcat(f{~cellfun(@isempty,f)});

% x = regexp(x{1,1},'[\.\-0-9]*','match');
v = cell2mat(textscan(v{1,1},'%*c %f %f %f'));
f = cell2mat(textscan(f{1,1},'%*c %f %f %f')); %it's really %d but