% function rotest(endname)
% clear
% endname = 'i2_540_t3_blob18k_a256e256L30';
endname = 'i2_540_t3_blob18k_a256e256L30r1e-06';
%% gather uniformly random (maybe uniform grid later) samples of unit quaternions i.e. S^3
if exist('SH/randq.mat','file')
  load SH/randq.mat
else
  N = 5e2;
%   N = 1e4;
  q = randn(N,4);
  q = bsxfun(@rdivide,q,sqrt(sum(q.^2,2)));
  save('SH/randq.mat','q','N');
end
%% load data to be rotationally matched
% load mcf/blob1k.mat
load mcf/blob18k.mat
vmcf = v;
load(['SH/' endname '.mat']);
%% compute cost of deviation at each random sample
J = zeros(N,1);
parfor i = 1:N
  J(i) = rotcost(q(i,:)',s_T,vmcf,s_end,v);
end
[J,idx] = sort(J);
q = q(idx,:);
%% locally optimize nop number of best candidates
nop = 5;
qqop = zeros(nop,4);
JJop = zeros(nop,1);
for i = 1:nop
  qq = q(i,:)';
  test = @(q) rotcost(q,s_T,vmcf,s_end,v);
  options = optimset('display','iter','maxiter',10000,'tolFun',1e-5,'tolx',1e-5);
  [qop,JJop(i)] = fminsearch(test,qq,options);
  qqop(i,:) = qop/norm(qop);
end
%% display result
% close all;
[Jop,i] = min(JJop);
% qinv = @(q) [q(1); -q(2:4)];
qop = qqop(i,:)';
qop = qop/norm(qop);
% qop = qinv(qop); % watch out for quat2dcm conventions / SU(2) double cover of SO(3)
v = quatrot(v',qop)';
v_end = quatrot(v_end',qop)';
%%
figh = visualizeSH(v,v_T,v_end,vmcf,f,f_T,[],s_end,s_T,...
J_hist,Jc_hist,D_0,D_T,D_endp,D_end,0);
%% rewrite history
hgexport(figh,['SH/' endname '.png'],...
  hgexport('factorystyle'), 'Format', 'png'); 
save(['SH/' endname '.mat'],'v','v_T','v_end','f','f_T','a_end','s_end','s_T',...
  'D_0','D_T','D_endp','D_end','J_hist','Jc_hist');
% end