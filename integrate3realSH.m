function realout = integrate3realSH(LM)
coeff = realSH2cpxSH(LM);
l = LM(:,1);
m = LM(:,2);
absm = abs(m);
ppp = absm.*[ 1; 1; 1];
ppn = absm.*[ 1; 1;-1];
pnp = absm.*[ 1;-1; 1];
pnn = absm.*[ 1;-1;-1];
npp = absm.*[-1; 1; 1];
npn = absm.*[-1; 1;-1];
nnp = absm.*[-1;-1; 1];
nnn = absm.*[-1;-1;-1];
realout = prod(coeff([1 2 3]))*int3cpxSH(l, ppp) + ...
          prod(coeff([1 2 6]))*int3cpxSH(l, ppn) + ...
          prod(coeff([1 5 3]))*int3cpxSH(l, pnp) + ...
          prod(coeff([1 5 6]))*int3cpxSH(l, pnn) + ...
          prod(coeff([4 2 3]))*int3cpxSH(l, npp) + ...
          prod(coeff([4 2 6]))*int3cpxSH(l, npn) + ...
          prod(coeff([4 5 3]))*int3cpxSH(l, nnp) + ...
          prod(coeff([4 5 6]))*int3cpxSH(l, nnn);
end
function coeff = realSH2cpxSH(LM)
m = LM(:,2);
coeff = zeros(size(LM));
for mi = 1:numel(m)
  if m(mi) == 0
    coeff(mi,:) = [1 0];
  elseif m(mi) < 0
    coeff(mi,:) = [-(-1)^m(mi)*1i/sqrt(2) 1i/sqrt(2)];
  else
    coeff(mi,:) = [(-1)^m(mi)/sqrt(2) 1/sqrt(2)];
  end
end
% coeff(m == 0,:) = repmat([1 0],numel(m(m == 0)),1);
% coeff(m > 0,:)  = [(-1).^m(m > 0)/sqrt(2)     repmat(1/sqrt(2), numel(m(m > 0)),1)];
% coeff(m < 0,:)  = [-(-1).^m(m < 0)*1i/sqrt(2) repmat(1i/sqrt(2),numel(m(m < 0)),1)];
  
end
function cpxout = int3cpxSH(l,m)
cpxout = sqrt(prod(2*l+1)/4/pi)*Wigner3j(l,m)*Wigner3j(l,[0;0;0]);
end