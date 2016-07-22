function [vham,fham] = xyz2hammer(v)
% turn cartesian coord in R^3 to its (triangulated) hammer projection \in [-1,1]^2

rv = sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
latv = acos(v(:,3)./rv) - pi/2;
lonv = atan2(v(:,2),v(:,1));
hamz = sqrt(1 + cos(latv).*cos(lonv/2));
hamx = cos(latv).*sin(lonv/2)./hamz;
hamy = sin(latv)./hamz;

fham = delaunayTriangulation([hamx hamy]);
fham = fham.ConnectivityList;

vham = [hamx hamy zeros(size(hamx))];
end