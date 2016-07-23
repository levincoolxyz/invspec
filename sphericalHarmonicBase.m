function realY = sphericalHarmonicBase(v,L)
% function realY = sphericalHarmonicBase(v,L)
% real spherical harmonic values on vertices of spherical mesh (v,f)
% 
% INPUT
% v     - [numv x 3] vertex positions of a spherical mesh
% L     - highest degree of spherical harmonics used
% 
% OUTPUT
% realY - [numv x numSH] real sphereical harmonic function values on vertices v

[AZ,EL,R] = cart2sph(v(:,1),v(:,2),v(:,3));

if std(R) >= 1e-2
  error('input mesh is not sufficiently sphere like ...');
end

realY = zeros(size(v,1),(L+1)^2);
SHidx = 0; % lazy linear index for spherical harmonics
for l = 0:L % degree
  for m = -l:l % order
    SHidx = SHidx+1;
    
    if m > 0
      part = @(x) sqrt(2)*(-1)^m*real(x);
    elseif m < 0
      part = @(x) sqrt(2)*(-1)^m*imag(x);
    else
      part = @(x) x;
    end
    m = abs(m);

    Pml = legendre(l,cos(EL + pi/2)); % Associated Legendre Polynomials of degree l
    if l~=0, Pml = Pml(m + 1,:); end % pick terms of order m

    C = sqrt((2*l + 1)/(4*pi) * gamma(l - m + 1)/gamma(l + m + 1)); % normalization constant
    Yml = C*Pml'.*exp(1i*m*AZ); % complex spherical harmonics
    realY(:,SHidx) = part(Yml); % convert to real SH basis
  end
end
end