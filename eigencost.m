function varargout = eigencost(s,M,L,lambda_T,numeig)
nums = numel(s);
%% get eigenvalues and eigenvectors
[V,lambda,ix] = eigvf(L,diag(1./s)*M,numeig);
%% (inversely weighted) mean square difference
% lambda_diff = lambda-lambda_T;
lambda_diff = 1./lambda-1./lambda_T;
% lambda_diff = 1./(lambda-lambda_T);
J = .5*sum(lambda_diff(1:(end-1)).^2); 
if (nargout <= 1)
  varargout{1} = J;
else
  varargout{1} = J;
  GJ = zeros(nums,1);
  for j = 1:nums
    for i = 1:(numeig-1)
      vi = V(:,ix(i));
      wi = zeros(nums,1);
%       wi(j) = lambda_diff(i)*lambda(i)*vi(j)^2/s(j);
      wi(j) = -lambda_diff(i)/lambda(i)*vi(j)^2/s(j);
%       wi(j) = lambda_diff(i)^(-3)*lambda(i)*vi(j)^2/s(j);
      GJ = GJ + wi;
    end
  varargout{2} = GJ;
  end
end