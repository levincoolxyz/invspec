function varargout = eigencostSH(a,mu_T,numeig)
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

% numa = numel(a);
numa = numeig;
if numeig > numel(mu_T), numeig = numel(mu_T); end
%% get eigenvalues and eigenvectors
[V,mu] = eigvfSH(a,numeig);
% [V,mu,L] = eigvfSH(a,numeig);
% Lstruct = L ./ L';
%% spectral cost and gradient
mu_T = mu_T(end-numeig+1:end); % chop unused target values
mu_diff = mu-mu_T;             % take the normal difference
mu_diff = mu_diff./mu_T;       % change to relative difference
% mu_diff = 1./mu-1./mu_T;       % take the inverse difference

J = .5*sum(mu_diff(1:(end-1)).^2); % compute norm squared cost

if (nargout <= 1)
  varargout{1} = J;
else
  varargout{1} = J;
  GJ = zeros(numa,1);
%   for i = 1:(numeig-1) % iterate relevant eigenvalue
%     vi = V(:,i);
%     mdi = mu_diff(i);
%     mTi = mu_T(i);
% %     mTi = -mu(i)^2; % inverse difference
%     parfor j = 1:numa % iterate relevant conformal factor
%       x = load(num2str(j,'../RSHIs/dRSHI%04d.mat'));
%       wij = mdi*vi'*x.pLpaj(1:numeig,1:numeig)*vi;   % normal difference
%       wij = wij/mTi;                                 % relative difference
%       GJ(j) = GJ(j) + wij;
%     end
%   end
  parfor j = 1:numa
    x = load(num2str(j,'../RSHIs/dRSHI%04d.mat'));
    dmuidaj = zeros(numeig-1,1);
    for i = 1:(numeig-1)
      dmuidaj(i) = V(:,i)'*x.pLpaj(1:numeig,1:numeig)*V(:,i);
    end
    GJ(j) = sum(mu_diff(1:end-1)./mu_T(1:end-1).*dmuidaj);
  end
  varargout{2} = GJ;
end