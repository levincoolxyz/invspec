function [a,xsh] = sphericalHarmonicDecomposition(v,x,L)
% function [a,xsh] = sphericalHarmonicDecomposition(v,f,x,n)
% spherical harmonic decomposition of scalar field x on spherical mesh (v,f)
% 
% INPUT
% v   - [numv x 3] vertex positions of a spherical mesh
% x   - [numv x 1] per-vertex scalar field on mesh (v,f)
% L   - highest degree of spherical harmonics used
% 
% OUTPUT
% a   - coefficient of the spherical harmonic decompostion, i.e.
%       x ~= sum(a(i)*real(Ymn).^2)
% xsh - the spherical harmonic approximated scalar field

vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);

[AZ,EL,R] = cart2sph(v(:,1),v(:,2),v(:,3));

if std(R) >= 1e-2
  error('input mesh is not sphere like (enough)');
end

% get (real) spherical harmonic functions
parfor vi = 1:size(v,1)
  idx = 0;
  az = AZ(vi);
  el = EL(vi) + pi/2;
  for l = 0:L
    for m = 0:l
      idx = idx+1;
      Pml = legendre(l,cos(el));

      if l~=0
        Pml = squeeze(Pml(m + 1,:,:));
      end

      a1 = (2*l + 1)/(4*pi);
      a2 = gamma(l - m + 1)/gamma(l + m + 1);
      C = sqrt(a1*a2);

      Yml = C*Pml.*exp(1i*m*az);

      realY{vi}(idx) = real(Yml);
      lambda{vi}(idx) = -l*(l+1);
    end
  end
end

realY = cell2mat(realY');
lambda = cell2mat(lambda');

% a = (realY'*realY)\(realY'*x);
a = lambda'*x;
xsh = realY*a;
end