clear;
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

% endname = 'i2_540_t3_blob1k_e49SH';
% endname = 'i4_mcf_t3_blob1k_e49SH';

% endname = 'i2_540_t3_blob18k_e36SH';
% endname = 'i2_540_t3_blob18k_e64SH';
% endname = 'i2_540_t3_blob18k_e49PL';
% endname = 'i2_540_t3_blob18k_e49L8';
% endname = 'i2_540_t3_blob18k_e49SH';
% endname = 'i2_540_t3_blob18k_e49exp';
% endname = 'i2_540_t3_blob18k_e49inv';
% endname = 'i2_540_t3_blob18k_a49e49L30';
endname = 'i2_540_t3_blob18k_a49e49L30s';
% endname = 'i2_540_t3_blob18k_a49e49L12';
% endname = 'i2_540_t3_blob18k_a64e64L30';
% endname = 'i2_540_t3_blob18k_a36e36L30';
% endname = 'i2_1000_t3_blob18k_e49SH';
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
close all;
[Jop,i] = min(JJop);
% qinv = @(q) [q(1); -q(2:4)];
qop = qqop(i,:)';
qop = qop/norm(qop);
% qop = qinv(qop); % watch out for quat2dcm conventions / SU(2) double cover of SO(3)
v = quatrot(v',qop)';
v_end = quatrot(v_end',qop)';
%%
figh = visualizeSH(v,v_T,v_end,vmcf,f,f_T,s_end,s_T,...
J_hist,Jc_hist,D_0,D_T,D_endp,D_end,0);
%% rewrite history
hgexport(figh,['SH/' endname '.png'],...
  hgexport('factorystyle'), 'Format', 'png'); 
save(['SH/' endname '.mat'],'v','v_T','v_end','f','f_T','s_end','s_T',...
  'D_0','D_T','D_endp','D_end','J_hist','Jc_hist');