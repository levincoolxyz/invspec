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
if numeig > numel(lambda_T), numeig = numel(lambda_T); end
%% get eigenvalues and eigenvectors
% [V,lambda] = eigvf(L,diag(1./s)*M,numeig);
[V,lambda] = eigvf(L,diag(1./exp(s))*M,numeig); % log-conformal factors
%% spectral cost and gradient
lambda_T = lambda_T(end-numeig+1:end); % chop unused target values
lambda_diff = lambda-lambda_T; % take the normal difference
lambda_diff = lambda_diff./lambda_T; % change to relative difference
% lambda_diff = 1./lambda-1./lambda_T; % take inverse difference

J = .5*sum(lambda_diff(1:(end-1)).^2); % compute norm squared cost

LLs = L*L*s;
J = J + s'*LLs/2*reg;  % apply bi-laplacian regularization
if (nargout <= 1)
  varargout{1} = J;
else
  varargout{1} = J;
  GJ = zeros(nums,1);
  m = diag(M);
  parfor j = 1:nums % iterate relevant conformal factor
    sj = s(j);
    Mjj = m(j);
    for i = 1:(numeig-eig0+1) % iterate relevant eigenvalue
      vij = V(j,i);
%       wij = lambda_diff(i)*lambda(i)*vij^2/sj^2*Mjj;   % normal difference
      wij = lambda_diff(i)*lambda(i)*vij^2*exp(-sj)*Mjj; % normal difference (log)
      wij = wij/lambda_T(i);                             % relative difference
%       wij = -lambda_diff(i)/lambda(i)*vij^2/sj^2*Mjj;  % inverse difference
%       wij = -lambda_diff(i)/lambda(i)*vij^2*exp(-sj)*Mjj; % inv difference (log)
      GJ(j) = GJ(j) + wij;
    end
    GJ(j) = GJ(j) + LLs(j)*reg;
  end
  varargout{2} = GJ;
end