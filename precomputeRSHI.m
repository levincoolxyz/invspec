clear all;
maxL = 100;
numeig = (maxL+1)^2;
%% compute full triple Real Spherical Harmonics Integral

[~,LM] = sphericalHarmonicBase([1 1 1],maxL);
for i = 1:numeig
  cijk = zeros(numeig);
  for j = 1:numeig
    akcijk = 0;
    parfor k = 1:numeig
      lm = LM([i j k],:);
      if ( lm(3) > (lm(1) + lm(2)) ) || ( lm(3) < abs(lm(1) - lm(2)) )
        cijk(j,k) = 0;
      else
        cijk(j,k) = integrate3realSH(lm);
      end
    end
  end
  cijk = real(cijk);
  save(num2str(i,'../RSHI/RSHI%04d.mat'),'cijk'); % need to buy more RAM
end
% save('../RSHI/RSHI.mat','cijk');

%% sparsify
for i = 1:numeig
  load(num2str(i,'../RSHI/RSHI%04d.mat'));
  cijk = sparse(cijk);
  save(num2str(i,'../RSHIs/RSHI%04d.mat'),'cijk'); % need to buy more RAM
end

%% efficient storage scheme?

%% visualize sparsity
for i = 1:numeig
  load(num2str(i,'../RSHI/RSHI%04d.mat'));
  spy(cijk)
  pause(.1)
end

%% precompute derivatives
% numeig = 36; %testing
D_s = [];
for l = 0:sqrt(numeig)-1
  D_s = [D_s, repmat(-l*(l+1),1,2*l+1)]; % spherical harmonic eigenvalues
end

% a = sym('a%d', [numeig 1]);
% assume('a%d', 'real');
% L = cell(numeig,1);
% for i = 1:numeig
%   load(num2str(i,'../RSHI/RSHI%04d.mat'));
%   L{i} = D_s(i)*a'*cijk(1:numeig,1:numeig);
% end

% pLpajcell = cell(numeig,1);
pLpaj = zeros(numeig);
for j = 11:numeig
  parfor i = 1:numeig
%     pLpajcell{i} = diff(L{i},a(j));
%     pLpaj(i,:) = double(pLpajcell{i});
    x=load(num2str(i,'../RSHIs/RSHI%04d.mat'));
%     da = sparse(1,j,1,1,numeig);
%     pLpaj(i,:) = D_s(i)*da*cijk;
%     pLpaj(i,:) = D_s(i)*da*cijk(1:numeig,1:numeig);
    pLpaj(i,:) = D_s(i)*x.cijk(j,:);
%     pLpaj(i,:) = D_s(i)*cijk(j,1:numeig);
  end
%   pLpaj = sparse(pLpaj);
  save(num2str(j,'../RSHIs/dRSHI%04d.mat'),'pLpaj'); % need to buy more RAM
end