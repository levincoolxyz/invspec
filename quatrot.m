function [w] = quatrot(v,q)
% function [w] = quatrot(v,q)
% quaternion rotation of (0+v) by q
% 
% INPUTS
% v       - [3 x n] array of vectors in \R^3
% s       - [1 x 1] cos( angle of rotation /2)
% t       - [3 x 1] unit vector axis of rotation in \R^3 * sin( angle of rotation /2)
% 
% OUTPUT
% w       - [3 x n] array of rotated vectors in \R^3

n = size(v,2);
q = repmat(q,1,n);
[u,c] = qmult(q(1,:),q(2:end,:),zeros(1,n),v);
[w,~] = qmult(c,u,q(1,:),-q(2:end,:));
end

function [u,c] = qmult(a,v,b,w)
% function [u,c] = qmult(a,v,b,w)
% quaternion multiplication of (a+v) and (b+w) to (c+u)

c = a.*b - dot(v,w,1);
u = cross(v,w,1) + repmat(a,3,1).*w + repmat(b,3,1).*v;
end

% function [w] = quatrot(v,u,theta)
% % function [w] = quatrot(v,u,theta)
% % quaternion rotation of v by angle theta around (unit) vector u
% % 
% % INPUTS
% % v       - [3 x n] array of vectors in \R^3
% % u       - [3 x 1] unit vector in \R^3
% % theta   - [1 x 1] angle of rotation (radian)
% % 
% % OUTPUT
% % w       - [3 x n] array of rotated vectors in \R^3
% 
% n = size(v,2);
% s = cos(theta/2);
% t = repmat(u*sin(theta/2),1,n);
% w = qmult(0,qmult(s,t,0,v),s,-t);
% end