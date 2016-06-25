function varargout = eigencost(s,M,L,lambda_T,numeig,reg,eig0)
% function [J,GJ] = eigencost(s,M,L,lambda_T,numeig,reg,eig0)
% 
% INPUT
% s        - (numv x 1) vector of per-vertex conformal factors
% M, L     - (numv x numv) mass and cotan matrices
% lambda_T - (numv x 1) target discrete spectrum
% numeig   - number of eigenvalues considered
% reg      - coefficient for the regularization term
% eig0     - the starting eigenvalue being considered
% OUTPUT
% J        - cost
% GJ       - gradient
% 

if nargin<6 || isempty(reg), reg = 0; end
if nargin<7 || isempty(eig0) || eig0<2, eig0 = 2; end
nums = numel(s);
%% get eigenvalues and eigenvectorsh
% [V,lambda] = eigvf(L,diag(1./s)*M,numeig);
[V,lambda] = eigvf(L,diag(1./exp(s))*M,numeig); % log-conformal factors
%% (weighted) difference squared as costs
lambda_T = lambda_T(end-numeig+1:end);
lambda_diff = lambda-lambda_T;
lambda_diff = lambda_diff./lambda_T;
% lambda_diff = 1./lambda-1./lambda_T;
J = .5*sum(lambda_diff(1:(end-1)).^2);
% Ls = L*s;
% J = J - s'*Ls/2*reg;
LLs = L*L*s;
J = J + s'*LLs/2*reg;
if (nargout <= 1)
  varargout{1} = J;
else
  varargout{1} = J;
  GJ = zeros(nums,1);
  parfor j = 1:nums
    for i = 1:(numeig-eig0+1)
      vi = V(:,i);
%       wij = lambda_diff(i)*lambda(i)*vi(j)^2/s(j)^2*M(j,j);
      wij = lambda_diff(i)*lambda(i)*vi(j)^2*exp(-s(j))*M(j,j);
      wij = wij/lambda_T(i);
%       wij = -lambda_diff(i)/lambda(i)*vi(j)^2/s(j)^2*M(j,j);
%       wij = -lambda_diff(i)/lambda(i)*vi(j)^2*exp(-s(j))*M(j,j);
      GJ(j) = GJ(j) + wij;
    end
%       GJ(j) = GJ(j) - Ls(j)*reg;
      GJ(j) = GJ(j) + LLs(j)*reg;
  end
  varargout{2} = GJ;
end