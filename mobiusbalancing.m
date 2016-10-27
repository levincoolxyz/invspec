function [vbalanced,c] = mobiusbalancing(v,v_T,f)
  numf = size(f,1);
  
%   totalInitialArea = totalArea(v_T,f);
%   initialArea = getFaceAreas(v_T,f);

  %initial face areas
%   tildeAijk = initialArea; 
  tildeAijk = getFaceAreas(v_T,f);
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
  c = X\b;
  vbalanced = sphinv(v,c');
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