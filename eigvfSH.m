function varargout = eigvfSH(a,numeig,maxL)
% [D] = eigvf(a,numeig,maxL)
% [V,D] = eigvf(a,numeig,maxL)
% eigenvalue/vector calculation in Spherical Harmonic Basis
% results sorted by ascending eigenvalue magnitude
% use numeig<=256 for 32 digits vpa

% if mod(sqrt(numeig),1) ~= 0
%   error('numeig must be a square for symmetry preservation');
% end

numL = (maxL+1)^2;
if numel(a)<numL
  a = [a;zeros(numL-numel(a),1)];
end

D_s = [];
for l = 0:sqrt(numL)-1
    D_s = [D_s, repmat(-l*(l+1),1,2*l+1)]; % spherical harmonic eigenvalues
end

L = zeros(numL);
for i = 1:numL
  if exist(num2str(i,'../RSHIs/RSHI%04d.mat'),'file')
    load(num2str(i,'../RSHIs/RSHI%04d.mat'));
  else
    error('need to precompute the integrals first');
  end
  L(i,:) = D_s(i)*a(1:numL)'*cijk(1:numL,1:numL);
end

if (nargout <= 1)
%   varargout{1} = sort(eig(L));
  lambda = sort(eig(L));
  varargout{1} = lambda(end-numeig+1:end);
else
  [V,D,W] = eig(L);
%   [varargout{2},I] = sort(diag(D));
%   varargout{1} = V(:,I);
%   varargout{3} = W(:,I);
  [lambda,I] = sort(real(diag(D)));
  varargout{2} = (lambda(end-numeig+1:end));
%   varargout{1} = V(end-numeig+1:end,I);
%   varargout{3} = W(end-numeig+1:end,I);
  varargout{1} = real(V(:,I(1:numeig)));
  varargout{3} = real(W(:,I(1:numeig)));
end