function realY = sphericalHarmonicBase(v,L)
% function realY = sphericalHarmonicBase(v,L)
% real spherical harmonic values on vertices of spherical mesh (v,f)
% 
% INPUT
% v     - [numv x 3] vertex positions of a spherical mesh
% L     - highest degree of spherical harmonics used
% 
% OUTPUT
% realY - [numv x numSH] (real) sphereical harmonic function values on vertices of v

[AZ,EL,R] = cart2sph(v(:,1),v(:,2),v(:,3));

if std(R) >= 1e-2
  error('input mesh is not sufficiently sphere like ...');
end

% get (real) spherical harmonic functions
parfor vi = 1:size(v,1)
  idx = 0;
  az = AZ(vi);
  el = EL(vi) + pi/2;
  for l = 0:L
    for m = -l:l
      if m > 0
        part = @(x) sqrt(2)*(-1)^m*real(x);
      elseif m < 0
        part = @(x) sqrt(2)*(-1)^m*imag(x);
      else
        part = @(x) x;
      end
      
      m = abs(m);
      
      idx = idx+1;
      Pml = legendre(l,cos(el));

      if l~=0
        Pml = squeeze(Pml(m + 1,:,:));
      end

      a1 = (2*l + 1)/(4*pi);
      a2 = gamma(l - m + 1)/gamma(l + m + 1);
      C = sqrt(a1*a2);

      Yml = C*Pml.*exp(1i*m*az); % complex spherical harmonics

      realY{vi}(idx) = part(Yml); % convert to real SH basis
    end
  end
end

realY = cell2mat(realY');
end