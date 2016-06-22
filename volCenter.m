function [zentrum,vol] = volCenter(v,f)
zentrum = zeros(1,3);
vol = 0;
for fi = 1:size(f,1)
  v1 = v(f(fi,1),:);
  v2 = v(f(fi,2),:);
  v3 = v(f(fi,3),:);
  vol = vol + det([v1;v2;v3]);
  zentrum = zentrum + (v1 + v2 + v3)*det([v1;v2;v3])/4;
end
zentrum = zentrum/vol;
end