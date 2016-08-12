clear all

maxL = 30;
[~,LM] = sphericalHarmonicBase([1 1 1],maxL);
numeig = size(LM,1);

for i = 117:numeig
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