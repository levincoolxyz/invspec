clear all

maxL = 30;
[~,LM] = sphericalHarmonicBase([1 1 1],maxL);
numeig = size(LM,1);
%%

for i = 856:numeig
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

%% memory test
CIJK = cell(numeig,1);
for i = 1:numeig
  load(num2str(i,'../RSHI/RSHI%04d.mat'));
  CIJK{i} = cijk;
end

%% visualize sparsity
for i = 1:numeig
  load(num2str(i,'../RSHI/RSHI%04d.mat'));
  spy(cijk)
  pause(.1)
end