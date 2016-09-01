function varargout = eigvfSH(a,numeig)
% [D] = eigvf(a,numeig)
% [V,D] = eigvf(a,numeig)
% eigenvalue/vector calculation in Spherical Harmonic Basis
% results sorted by ascending eigenvalue magnitude
% use numeig<=256 for 32 digits vpa

if mod(sqrt(numeig),1) ~= 0
  error('numeig must be a square for symmetry preservation');
end

D_s = [];
for l = 0:sqrt(numeig)-1
    D_s = [D_s, repmat(-l*(l+1),1,2*l+1)]; % spherical harmonic eigenvalues
end

L = zeros(numeig);
for i = 1:numeig
  if exist(num2str(i,'../RSHI/RSHI%04d.mat'),'file')
    load(num2str(i,'../RSHI/RSHI%04d.mat'));
  else
    error('need to precompute the integrals first');
  end
  L(i,:) = D_s(i)*a'*cijk(1:numeig,1:numeig);
end

if (nargout <= 1)
  varargout{1} = sort(eig(L));
else
  [V,D] = eig(L);
  [varargout{2},I] = sort(D);
  varargout{1} = V(:,I);
end