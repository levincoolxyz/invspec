function varargout = eigencost(s,M,L,lambda_T,numeig)
% function [J,GJ] = eigencost(s,M,L,lambda_T,numeig)
nums = numel(s);
%% get eigenvalues and eigenvectors
[V,lambda] = eigvf(L,diag(1./s)*M,numeig);
%% (weighted) difference squared as costs
% lambda_diff = lambda-lambda_T;
lambda_diff = 1./lambda-1./lambda_T;
J = .5*sum(lambda_diff(1:(end-1)).^2); 
if (nargout <= 1)
  varargout{1} = J;
else
  varargout{1} = J;
  GJ = zeros(nums,1);
  parfor j = 1:nums
    for i = 1:(numeig-1)
      vi = V(:,i);
%       wij = lambda_diff(i)*lambda(i)*vi(j)^2/s(j)^2*M(j,j);
      wij = -lambda_diff(i)/lambda(i)*vi(j)^2/s(j)^2*M(j,j);
      GJ(j) = GJ(j) + wij;
    end
  end
  varargout{2} = GJ;
end