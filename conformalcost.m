function [J,GJ] = conformalcost(v,isedge,elsq_T)
%% initialization
v = reshape(v,3,[])';
numv = size(v,1); % number of vertices
%% compute (upper triangular) edge lengths and squared
el = zeros(numv,3);
elsq = zeros(numv);
for i = 1:numv
%   for j = find(isedge(i,:))
  for j = (i+1):numv % linear indices
    for dim = 1:3
      el(i,j,dim) = v(i,dim)-v(j,dim);
    end
    elsq(i,j) = sum((v(i,:)-v(j,:)).^2);
  end
end
elsq = elsq(isedge); % linear indices
%% naively compare edge lengths squared difference, cost, and gradient
elsq_diff = elsq - elsq_T;
% J = .25*sum(sum(elsq_diff.^2));
J = .25*sum(elsq_diff.^2);
GJ = zeros(numv,3);
for vi = 1:numv
  for dim = 1:3
    ddvi = zeros(numv);
    ddvi(vi,:) = el(vi,:,dim);
    ddvi(:,vi) = -el(:,vi,dim);
%     GJ(vi,dim) = sum(sum(elsq_diff.*ddvi));
    GJ(vi,dim) = sum(elsq_diff.*ddvi(isedge));
  end
end
GJ = reshape(GJ',[],1);
%% using spring analogy, cost, and gradient
% elsq_rat = elsq_T./elsq;
% J = .25*sum((1 - elsq_rat).^2);
% GJ = zeros(numv,3);
% for vi = 1:numv
%   for dim = 1:3
%     ddvi = zeros(numv);
%     ddvi(vi,:) = el(vi,:,dim);
%     ddvi(:,vi) = -el(:,vi,dim);
%     GJ(vi,dim) = sum((1 - elsq_rat).*elsq_T./elsq./elsq.*ddvi(isedge));
%   end
% end
% GJ = reshape(GJ',[],1);