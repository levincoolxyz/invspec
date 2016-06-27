function varargout = eigvf(L,M,numeig)
% [D] = eigvf(L,M,numeig)
% [V,D] = eigvf(L,M,numeig)
% custom generalized eigenvalue problem wrapper function
% results sorted by ascending eigenvalue magnitude

if (nargout <= 1)
  if numeig < size(M,1)
    D = eigs(L,M,numeig,-1e-6);
  elseif ~issparse(L)
    D = eig(L,M);
  else
    warning('ignoring highest frequency; matlab is a bitch');
    D = eigs(L,sparse(M),numeig-1,-1e-6);
  end
  varargout{1} = sort(D);
else
  norm = @(v) sqrt(dot(v,M*v));
  if numeig < size(M,1)
    [V,D] = eigs(L,M,numeig,-1e-6);
    % normalize eigenvectors
    for i = 1:numeig
      V(:,i) = V(:,i)./norm(V(:,i));
    end
  elseif ~issparse(L)
    [V,D] = eig(L,M);
    % normalize eigenvectors
    for i = 1:size(M,1)
      V(:,i) = V(:,i)./norm(V(:,i));
    end
  else
    warning('ignoring highest frequency; matlab is a bitch');
    D = eigs(L,sparse(M),numeig-1,-1e-6);
  end
  [varargout{2},I] = sort(diag(D));
  varargout{1} = V(:,I);
end