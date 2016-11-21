function [vbalanced,c] = mobiusbalancing(v,v_T,f)
  numf = size(f,1);
  numv = size(v,1);
  if numv ~= size(v_T,1), error('vertex sets size mismatch'); end
  v0 = v;

  %target face areas
  tildeAijk = getFaceAreas(v_T,f);
%   totalInitialArea = totalArea(v_T,f);

  vnew = zeros(size(v));
%   imax = 1e4;
  imax = 2e3;
  iter = 1;
  c = [0;0;0];
  while iter < imax
    %current face areas
    Aijk = getFaceAreas(v,f);
    % current face barycenters
    fijk = getFaceBarycenters(v,f);

    % current surface area
    A = totalArea(v,f);

    % Compute the square of the conformal factor per face 
    phi = tildeAijk./Aijk/2;

    % Compute the mean log scale factor 
    u0 = sum(Aijk.*log(sqrt(phi)))/A;

    % Compute target flow direction xi for squared conformal factors 
    phi0 = exp(2*u0);
    xi = phi - phi0;

    % Compute right-hand side of Euler-Lagrange system
    b = -sum(repmat(Aijk.*xi,1,3).*fijk,1)'/2;

    % Compute matrix for Euler-Lagrange system
    X = zeros(3,3);
    for j = 1:numf
      X = X + Aijk(j) * ( fijk(j,:)' * fijk(j,:));
    end
    a = (X\b);
    c = c + a;
    vnew = sphinv(v0,c);
%     scl = makeUnitArea(vnew,f);
%     vnew = vnew - repmat(volCenter(vnew,f),numv,1);
%     vnew = vnew*scl*sqrt(4*pi);
    iter = iter + 1;
%     mean(sum(abs(vnew - v)./abs(vnew),2))
%     norm(a)
    if mean(sum(abs(vnew - v)./abs(vnew),2)) <= .1
      break;
    end
    tildeAijk = getFaceAreas(v,f);
    v = vnew;
%     if iter > 500
%     trimesh(f,v(:,1),v(:,2),v(:,3));
%     end
  end
  
  vbalanced = vnew;
end
function area = getFaceAreas(v,f)
  area = zeros(size(f,1),1);
  for fi = 1:size(f,1)
    vi = f(fi,:);
    area(fi) = norm(cross(v(vi(2),:) - v(vi(1),:), v(vi(3),:) - v(vi(1),:)))/2;
  end
end
function bary = getFaceBarycenters(v,f)
bary = zeros(size(f));
for fi = 1:size(f,1)
  vi = f(fi,:);
  bary(fi,:) = sum(v(vi,:),1)/3;
end
end
function totalarea = totalArea(v,f)
  totalarea = sum(getFaceAreas(v,f));
end