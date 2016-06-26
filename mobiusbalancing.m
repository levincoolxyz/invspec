function v_F = mobiusbalancing(v,f,s0,s_T)
  totalInitialArea = makeUnitArea(v,f);

  %initial face areas
  for j = 1:numf
    tildeAijk(j) = initialFaceArea;
  %current face areas
    Aijk(j) = getFaceArea(v,f);
  % current face barycenters
  fijk = getFaceBarycenters(v,f);
  end
  
  % current surface area
  A = makeUnitArea(v,j);

  % Compute the square of the conformal factor per face 
  phi = zeros(numf,1);
  for j = 1:numf
    phi(j) = (1./2.) .* tildeAijk(j) ./ Aijk(j);
  end

  % Compute the mean log scale factor 
  u0 = 0;
  for j = 1:numf
    u0 = u0 + Aijk(j) .* log( sqrt( phi(j) ));
  end
  u0 = u0 / A;

  % Compute target flow direction xi for squared conformal factors 
  xi = zeros(numf,3);
  phi0 = exp( 2*u0 );
  for j = 1:numf
    xi(j) = phi(j) - phi0;
  end

  % Compute right-hand side of Euler-Lagrange system 
  b = zeros(3,1);
  for j = 1:numf
    b = b + Aijk(j) * xi(j) * fijk(j);
  end
  b = b*-1/2;

  % Compute matrix for Euler-Lagrange system
  X = zeros(3,3);
  for j = 1:numf
    X = X + Aijk(j) * ( fijk(:,j) * fijk(:,j)');
  end
  a = X\b;
  v_F = sphinv(v,a);
end