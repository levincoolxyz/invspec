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

%% initialize and compute edge lengths with squares
v = reshape(v,3,[])';
numv = size(v,1); % number of vertices
elsq = zeros(numv); % edge length squares
col = mod(isedge-1,numv)+1; % isedge column indices
row = fix((isedge-1)/numv); % isedge row indices
if (nargout <= 1)
  for j = 1:numv
    for i = col(row == (j-1))'
      elsq(i,j) = sum((v(i,:)-v(j,:)).^2);
    end
  end
  elsq = nonzeros(elsq);
else
  el = zeros(numv,numv,3);
  for j = 1:numv
    for i = col(row == (j-1))'
      for dim = 1:3
        el(i,j,dim) = v(i,dim)-v(j,dim);
      end
      elsq(i,j) = sum((v(i,:)-v(j,:)).^2);
    end
  end
  elsq = nonzeros(elsq);
end
%% edge lengths cost and gradient
elsq_diff = elsq - elsq_T; % normal difference
% elsq_diff = elsq_diff./elsq; % relative difference (technically more correct)
J = .25*sum(elsq_diff.^2); % compute norm squared cost

if (nargout <= 1)
  varargout{1} = J;
else
  varargout{1} = J;
  GJ = zeros(numv,3);
  elcpy = el;
  parfor vi = 1:numv
    for dim = 1:3
      ddvi = zeros(numv);
      ddvi(vi,:) = el(vi,:,dim);
      ddvi(:,vi) = -elcpy(:,vi,dim);
      GJ(vi,dim) = sum(elsq_diff.*ddvi(isedge));
%       GJ(vi,dim) = sum(elsq_diff.*elsq_T./elsq./elsq.*ddvi(isedge));
    end
  end
  varargout{2} = reshape(GJ',[],1);
end