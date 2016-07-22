function M = padcat(varargin)
% padcat - concatenate 2D matrices of different row length by padding with NaN

varargin = varargin(~cellfun('isempty',varargin)); 

for i = 1:numel(varargin)
  rsz(i) = size(varargin{i},1);
  csz(i) = size(varargin{i},2);
  if rsz(i) == 0 && csz(i) == 0, varargin{i} = []; end
end

maxrsz = max(rsz);

for i = 1:numel(varargin)
  varargin{i} = [varargin{i}; nan(maxrsz-rsz(i),csz(i))];
end

M = cell2mat(varargin);