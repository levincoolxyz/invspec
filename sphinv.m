function vinv = sphinv(v,c)
% spherical inversion of vertices v at center c
% 
% INPUT
% v    - numv x 3 matrix of original vertex positions
% c    - 1 x 3 row vector of inversion center
% OUTPUT
% vinv - numv x 3 matrix of inverted vertex positions
% 

vnorm2 = @(v) (v(:,3).^2+v(:,1).^2+v(:,2).^2);
numv = size(v,1);
vinv = v./repmat(vnorm2(v),1,3)-repmat(c,numv,1);
vinv = vinv./repmat(vnorm2(vinv),1,3);
end