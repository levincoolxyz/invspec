function varargout = eigvf(L,M,numeig)

if (nargout <= 1)
  if numeig < size(M,1)
    D = eigs(L,M,numeig,-1e-6);
  else
%     D = eig(inv(M)*L);
%     D = eig(M\L);
    D = eig(L,M);
  end
  varargout{1} = sort(D);
else
%   Mnorm = @(v) sqrt(dot(v,M*v)); %does not work?
  Mnorm = @(v) sqrt(dot(v,v));
  if numeig < size(M,1)
    [V,D] = eigs(L,M,numeig,-1e-6);
    % normalize eigenvectors
    for i = 1:size(V,2)
      V(:,i) = V(:,i)./Mnorm(V(:,i));
    end
  else
%     [V,D] = eig(inv(M)*L);
%     [V,D] = eig(M\L);
    [V,D] = eig(L,M);
    % normalize eigenvectors
    for i = 1:size(V,2)
      V(:,i) = V(:,i)./Mnorm(V(:,i));
    end
  end
  varargout{1} = V;
  [varargout{2},varargout{3}] = sort(diag(D));
end