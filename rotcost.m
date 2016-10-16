function [J,GJ] = rotcost(q,s_T,v_T,s,v)
q = q/norm(q);
vnew = quatrot(v',q)';
I = zeros(size(s));
for vi = 1:size(v,1)
%   vnew = quat2dcm(q')*v(vi,:)'; %quat2dcm is different convention
%   [~,I(vi)] = max(dot(repmat(vnew',size(v_T,1),1),v_T,2));
  [~,I(vi)] = max(dot(repmat(vnew(vi,:),size(v_T,1),1),v_T,2));
end

J = 1/2*sum((s_T(I) - s).^2);

GJ = [0;0;0;0];
end

% function [J,GJ] = rotcost(q,sa_T,v_T,f_T,s,v,f)
% 
% vnew = quatrot(v',q)';
% s_T = zeros(size(s));
% for vi = 1:size(v,1)
%   s_T(vi) = sa_T(findface(vnew(vi,:),v_T,f_T));
% end
% 
% J = 1/2*sum((s_T - s).^2);
% 
% GJ = [0;0;0;0];
% end

function [idx,bc] = findface(p,v,f)
idx = 0;
bc = zeros(size(f));
for fi = 1:size(f,1)
  vi = f(fi,:);
%   if ifptintri3(p,v(vi,:))
  bc(fi,:) = ifptintri3(p,v(vi,:));
%   if sum(bc <= 1) ==3 && sum(bc >= 0) == 3
  if sum(abs(bc(fi,:)-.5) <= .5+1e-2) == 3
    idx = fi;
    break
  end
end
end

function answer = ifptintri3(p,pt)
% answer = ifptintri3(p0,p1,p2,p3)
% true if p is inside the triangle bdd by pt(1,:),pt(2,:),pt(3,:) if projected to its plane

% find p0's barycentric coordinate
u = pt(2,:) - pt(1,:);
v = pt(3,:) - pt(1,:);
w = p - pt(1,:);
n = cross(u,v);

bc(3) = dot(cross(u,w),n)/dot(n,n);
bc(2) = dot(cross(w,v),n)/dot(n,n);
bc(1) = 1 - bc(2) - bc(3);

answer = bc;
% answer = (sum(bc <= 1 + bc >= 0) == 6);
end