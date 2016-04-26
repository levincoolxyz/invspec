function [dcm,MoI] = pa(v,f)
% function [dcm,MoI] = pa(v,f)
% calculate principal frame and associated inertia of a mesh
% 
% INPUTS
% v,f    - face-vertex data of the mesh
% OUTPUTS
% dcm    - principal frame / direction cosine matrix for rotation
% MoI    - principal moments of inertia i.e. I_xx <= I_yy <= I_zz
% 
% Reference: 
% http://www.mjoldfield.com/atelier/2011/03/tetra-moi.html

numv = size(v,1);
numf = size(f,1);
[zentrum,vol] = volCenter(v,f);
% shift to volume center (and assume that it's inside the mesh)
v = v - repmat(zentrum,numv,1);

I = zeros(3);
for i = 1:numf
  t1 = 0;
  t2 = 0;
  t3 = 0;
  for j = 1:3
    vj = v(f(i,j),:);
    t1 = t1 + vj'*vj;
    t2 = t2 + vj';
    t3 = t3 + vj;
  end
  I = I + t1 + t2*t3;
end
[dcm,MoI] = eig(I);
[MoI,idx] = sort(diag(MoI));
MoI = MoI * vol / 120; % assume unit density
dcm = dcm(:,idx);
end

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