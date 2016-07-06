function varargout = eigvf(L,M,numeig)
% [D] = eigvf(L,M,numeig)
% [V,D] = eigvf(L,M,numeig)
% custom generalized eigenvalue problem wrapper function
% results sorted by ascending eigenvalue magnitude

sigma = -1e-6;

if (nargout <= 1)
  if numeig < size(M,1)
    D = eigs(L,M,numeig,sigma);
  elseif ~issparse(L)
    D = eig(L,M);
  else
    warning('ignoring highest frequency; matlab is a bitch');
    D = [NaN; eigs(L,sparse(M),numeig-1,sigma)];
  end
  varargout{1} = sort(D);
else
  norm = @(v) sqrt(dot(v,M*v));
  if numeig < size(M,1)
    [V,D] = eigs(L,M,numeig,sigma);
    % normalize eigenvectors
    for i = 1:numeig
      V(:,i) = V(:,i)./norm(V(:,i));
    end
    D = diag(D);
  elseif ~issparse(L)
    [V,D] = eig(L,M);
    % normalize eigenvectors
    for i = 1:size(M,1)
      V(:,i) = V(:,i)./norm(V(:,i));
    end
    D = diag(D);
  else
    warning('ignoring highest frequency; matlab is a bitch');
    [V,D] = eigs(L,sparse(M),numeig-1,sigma);
    % normalize eigenvectors
    for i = 1:size(M,1)
      V(:,i) = V(:,i)./norm(V(:,i));
    end
    D = [NaN; diag(D)];
  end
  [varargout{2},I] = sort(D);
  varargout{1} = V(:,I);
end