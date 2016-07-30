clear all

maxL = 30;
[~,LM] = sphericalHarmonicBase([1 1 1],maxL);
numeig = size(LM,1);

cijk = zeros(numeig);
for i = 688:700
  for j = 1:numeig
    akcijk = 0;
    parfor k = 1:numeig
      cijk(j,k) = integrate3realSH(LM([i j k],:));
    end
  end
  save(num2str(i,'../RSHI/RSHI%04d.mat'),'cijk'); % need to buy more RAM
end

% save('../RSHI/RSHI.mat','cijk');