clear all

maxL = 10;
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
  cijk = real(cijk);
  save(num2str(i,'../RSHI/RSHI%04d.mat'),'cijk'); % need to buy more RAM
end

% save('../RSHI/RSHI.mat','cijk');