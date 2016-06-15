function varargout = eigencost(s,M,L,lambda_T,numeig)
% function [J,GJ] = eigencost(s,M,L,lambda_T,numeig)
nums = numel(s);
%% get eigenvalues and eigenvectors
% [V,lambda] = eigvf(L,diag(1./s)*M,numeig);
[V,lambda] = eigvf(L,diag(1./exp(s))*M,numeig); % log-conformal factors
%% (weighted) difference squared as costs
% lambda_diff = lambda-lambda_T;
lambda_diff = 1./lambda-1./lambda_T;
J = .5*sum(lambda_diff(1:(end-1)).^2);
% Ls = L*s;
% J = J - s'*Ls/2;
if (nargout <= 1)
  varargout{1} = J;
else
  varargout{1} = J;
  GJ = zeros(nums,1);
  parfor j = 1:nums
    for i = 1:(numeig-1)
      vi = V(:,i);
%       wij = lambda_diff(i)*lambda(i)*vi(j)^2/s(j)^2*M(j,j);
%       wij = -lambda_diff(i)/lambda(i)*vi(j)^2/s(j)^2*M(j,j);
      wij = -lambda_diff(i)/lambda(i)*vi(j)^2*exp(-s(j))*M(j,j);
      GJ(j) = GJ(j) + wij;
    end
%       GJ(j) = GJ(j) - Ls(j);
  end
  varargout{2} = GJ;
end