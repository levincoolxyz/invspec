function h = plot3v(M,lspec)
numv = size(M,1);
v0 = zeros(numv,1);
h = quiver3(v0,v0,v0,M(:,1),M(:,2),M(:,3),lspec);