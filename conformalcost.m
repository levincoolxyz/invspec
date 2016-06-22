function varargout = conformalcost(v,isedge,elsq_T)
% function [J,GJ] = conformalcost(v,isedge,elsq_T)
% 
% INPUT
% v      - (numv x 3) vertex position coordiantes
% isedge - (nume x 1) linear indices of connected edges in mesh
% elsq_T - (nume x 1) target edge lengths
% OUTPUT
% J      - cost
% GJ     - gradient
% 

%% initialization
v = reshape(v,3,[])';
numv = size(v,1); % number of vertices
%% compute (upper triangular) edge lengths and squared
if (nargout > 1)
  el = zeros(numv,numv,3);
end
elsq = zeros(numv);
temp = mod(isedge-1,numv)+1;
temp2 = fix((isedge-1)/numv);
for j = 1:numv
  for i = temp(temp2 == (j-1))'
    if (nargout > 1)
      for dim = 1:3
        el(i,j,dim) = v(i,dim)-v(j,dim);
      end
    end
    elsq(i,j) = sum((v(i,:)-v(j,:)).^2);
  end
end
elsq = nonzeros(elsq);
%% edge lengths squared difference, cost, and gradient
elsq_diff = elsq - elsq_T;
% elsq_rat = elsq_T./elsq; % technically correct version
J = .25*sum(elsq_diff.^2);
% J = .25*sum((1 - elsq_rat).^2);
if (nargout <= 1)
  varargout{1} = J;
else
  varargout{1} = J;
  GJ = zeros(numv,3);
  parfor vi = 1:numv
    for dim = 1:3
      ddvi = zeros(numv);
      ddvi(vi,:) = el(vi,:,dim);
      ddvi(:,vi) = -el(:,vi,dim);
      GJ(vi,dim) = sum(elsq_diff.*ddvi(isedge));
%       GJ(vi,dim) = sum((1 - elsq_rat).*elsq_T./elsq./elsq.*ddvi(isedge));
    end
  end
  varargout{2} = reshape(GJ',[],1);
end