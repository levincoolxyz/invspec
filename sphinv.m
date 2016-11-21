function vinv = sphinv(v,c)
% double spherical inversion of vertices v at origin and center c
% 
% INPUT
% v    - numv x 3 matrix of original vertex positions \in \R^3
% c    - 1 x 3 row vector of inversion center \in \R^3
% OUTPUT
% vinv - numv x 3 matrix of inverted vertex positions
% 

c = c(:)';
numv = size(v,1);
vnorm2 = @(v) (v(:,1).^2+v(:,2).^2+v(:,3).^2);
vinv = -v./repmat(vnorm2(v),1,3)-repmat(c,numv,1);
vinv = -vinv./repmat(vnorm2(vinv),1,3);
end