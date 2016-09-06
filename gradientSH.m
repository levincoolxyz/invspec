clear all;
%%
numeig = 4; %testing
numa = numeig;
load '/home/ultimate/invspec/mcode/SH/bunnySHspecL=30.mat'
mu_T = eigvfSH(a_pj,numa);

D_s = [];
for l = 0:sqrt(numeig)-1
D_s = [D_s, repmat(-l*(l+1),1,2*l+1)]; % spherical harmonic eigenvalues
end
%%
a = sym('a%d', [numeig 1]);
assume(a, 'real');
L = sym('l',numeig);
for i = 1:numeig
  load(num2str(i,'../RSHI/RSHI%04d.mat'));
  L(i,:) = D_s(i)*a'*cijk(1:numeig,1:numeig);
end
%%
[V,mu] = eig(L);
mu = diag(mu);
%%
mu_T = mu_T(end-numeig+1:end); % chop unused target values
mu_diff = mu-mu_T;             % take the normal difference
mu_diff = mu_diff./mu_T;       % change to relative difference

J = .5*sum(mu_diff(1:(end-1)).^2); % compute norm squared cost
%%
for j = 1:numa
  x = load(num2str(j,'../RSHIs/dRSHI%04d.mat'));
  dmuidaj = sym('dmuidaj%d',[numeig-1,1]);
  for i = 1:(numeig-1)
    dmuidaj(i) = V(:,i)'*x.pLpaj(1:numeig,1:numeig)*V(:,i);
  end
  GJ(j) = sum(mu_diff(1:end-1)./mu_T(1:end-1).*dmuidaj);
end
%%
for j = 1:numa
  nablaJ(j) = diff(J,a(j));
end