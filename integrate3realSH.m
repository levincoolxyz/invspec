function realout = integrate3realSH(LM)
l = LM(:,1);
m = LM(:,2);
absm = abs(m);
mm = zeros(3,8);
mm(:,1) = absm.*[ 1; 1; 1]; %ppp
mm(:,2) = absm.*[ 1; 1;-1]; %ppn
mm(:,3) = absm.*[ 1;-1; 1]; %pnp
mm(:,4) = absm.*[ 1;-1;-1]; %pnn
mm(:,5) = absm.*[-1; 1; 1]; %npp
mm(:,6) = absm.*[-1; 1;-1]; %npn
mm(:,7) = absm.*[-1;-1; 1]; %nnp
mm(:,8) = absm.*[-1;-1;-1]; %nnn
mchk = sum(mm,1) == 0;
realSH = zeros(8,1);
for i = 1:8
  if mchk(i)
    realSH(i) = int3cpxSH(l, mm(:,i));
  end
end
coeff = realSH2cpxSH(LM);
realcoeff = [prod(coeff([1 2 3])) ...
             prod(coeff([1 2 6])) ...
             prod(coeff([1 5 3])) ...
             prod(coeff([1 5 6])) ...
             prod(coeff([4 2 3])) ...
             prod(coeff([4 2 6])) ...
             prod(coeff([4 5 3])) ...
             prod(coeff([4 5 6]))];
realout = realcoeff*realSH;
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