[X,Y,Z] = sphere(15);
% sphar = (2*Z.^2-X.^2-Y.^2);
sphar = 1./(abs(Z)+1);
% sphar = exp(-sphar);
% sphar = (3*X.^2-Y.^2).*Y.*Z;
surf(X,Y,Z,sphar)
colorbar