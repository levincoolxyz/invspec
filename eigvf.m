function varargout = eigvf(L,M,numeig)

% Mnorm = @(v) sqrt(dot(v,M*v)); %why does this ruin things?
Mnorm = @(v) sqrt(dot(v,v));

if numeig < size(M,1)
  [V,D] = eigs(L,M,numeig,-1e-6);
  % normalize eigenvectors
  for i = 1:size(V,2)
    V(:,i) = V(:,i)./Mnorm(V(:,i));
  end
else
%   [V,D] = eig(inv(M)*L);
%   [V,D] = eig(M\L);
  [V,D] = eig(L,M);
  % normalize eigenvectors
  for i = 1:size(V,2)
    V(:,i) = V(:,i)./Mnorm(V(:,i));
  end
end

if (nargout <= 1)
  varargout{1} = sort(diag(D));
else
  varargout{1} = V;
  [varargout{2},varargout{3}] = sort(diag(D));
end