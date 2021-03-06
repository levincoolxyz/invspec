function varargout = eigencostSH(a,mu_T,numeig,maxL,reg)
% function [J,GJ] = eigencostSH(a,lambda_T,numeig,maxL,reg)
% In Spherical Harmonics basis
% 
% INPUT
% a        - ((maxL+1)^2 x 1) vector of 1 / conformal factors in SH basis
% lambda_T - (numv x 1) target discrete spectrum
% numeig   - number of eigenvalues considered
% maxL     - max degree of SH basis considered
% reg      - bilaplacin regularization weight
% OUTPUT
% J        - cost
% GJ       - gradient
% 

numa = numel(a);
numL = (maxL+1)^2;
if numeig > numel(mu_T), numeig = numel(mu_T); end
%% get eigenvalues and eigenvectors
[V,mu,W,L] = eigvfSH(a,numeig,maxL);
% Lstruct = L ./ L'; % non symmetric matrix of only real eigenvalues IRL?!
%% spectral cost and gradient
mu_T = mu_T(end-numeig+1:end); % chop unused target values
mu_diff = mu-mu_T;             % take the normal difference
mu_diff = mu_diff./mu_T;       % change to relative difference
% mu_diff = 1./mu-1./mu_T;       % take the inverse difference

J = .5*sum(mu_diff(1:(end-1)).^2); % compute norm squared cost
L = L(1:numa,1:numa); % assume (maxL+1)^2>numa (needed for precision anyway)
La = L*a;
J = J + reg*(La'*La); % bilaplacian regularization

if (nargout <= 1)
  varargout{1} = J;
else
  varargout{1} = J;
  GJ = zeros(numa,1);
  parfor j = 1:numa
    x = load(num2str(j,'../RSHIs/dRSHI%04d.mat'));
    dmuidaj = zeros(numeig-1,1);
    pLpaj = x.pLpaj(1:numL,1:numL);
    for i = 1:(numeig-1)
      dmuidaj(i) = W(:,i)'*pLpaj*V(:,i)/(W(:,i)'*V(:,i));
    end
    GJ(j) = sum(mu_diff(1:end-1)./mu_T(1:end-1).*dmuidaj);
    pLapaj = pLpaj(1:numa,1:numa)*a+L*sparse(j,1,1,numa,1);
    GJ(j) = GJ(j) + reg*(pLapaj'*La+La'*pLapaj);
  end
  varargout{2} = GJ;
end

% to do #1: work out arbitrary (inverse) weighting case
% to do #2: work out repeating eigenvalue cases