function varargout = eigencostSH(a,lambda_T,numeig)
% function [J,GJ] = eigencostSH(a,lambda_T,numeig)
% In Spherical Harmonics basis
% 
% INPUT
% a        - ((maxL+1)^2 x 1) vector of 1 / conformal factors in SH basis
% lambda_T - (numv x 1) target discrete spectrum
% numeig   - number of eigenvalues considered (must be a square)
% OUTPUT
% J        - cost
% GJ       - gradient
% 

numa = numel(a);
if numeig > numel(lambda_T), numeig = numel(lambda_T); end
%% get eigenvalues and eigenvectors
[V,lambda] = eigvfSH(a,numeig);
%% spectral cost and gradient
lambda_T = lambda_T(end-numeig+1:end); % chop unused target values
lambda_diff = lambda-lambda_T;         % take the normal difference
lambda_diff = lambda_diff./lambda_T;   % change to relative difference

J = .5*sum(lambda_diff(1:(end-1)).^2); % compute norm squared cost

if (nargout <= 1)
  varargout{1} = J;
else
  varargout{1} = J;
  GJ = zeros(numa,1);
  for j = 1:numa % iterate relevant conformal factor
    aj = a(j);
    for i = 1:(numeig-eig0+1) % iterate relevant eigenvalue
      vij = V(j,i);
      wij = -lambda_diff(i)*lambda(i)*vij^2/aj^2;   % normal difference
      wij = wij/lambda_T(i);                        % relative difference
      GJ(j) = GJ(j) + wij;
    end
  end
  varargout{2} = GJ;
end