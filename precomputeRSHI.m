clear all

maxL = 50;
[~,LM] = sphericalHarmonicBase([1 1 1],maxL);
numeig = size(LM,1);

cijk = zeros(numeig);
for i = 1:numeig
  for j = 1:numeig
    akcijk = 0;
    parfor k = 1:numeig
      cijk(j,k) = integrate3realSH(LM([i j k],:));
    end
  end
  save(num2str(i,'../RSHI/RSHI%02d.mat'),'cijk'); % need to buy more RAM
end

% save('../RSHI/RSHI.mat','cijk');