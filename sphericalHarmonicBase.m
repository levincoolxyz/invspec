function [realY,LM] = sphericalHarmonicBase(v,L)
% function realY = sphericalHarmonicBase(v,L)
% real spherical harmonic values on vertices of spherical mesh (v,f)
% 
% INPUT
% v     - [numv x 3] vertex positions of a spherical mesh
% L     - highest degree of spherical harmonics used
% 
% OUTPUT
% realY - [numv x numSH] real sphereical harmonic function values on vertices v
% LM    - [numSH x 2] sphereical harmonic degrees and orders

[AZ,EL,R] = cart2sph(v(:,1),v(:,2),v(:,3));

if std(R) >= 1e-2
  error('input mesh is not sufficiently sphere like ...');
end

LM = [];
numv = size(v,1);
realY = zeros(numv,(L+1)^2);
for l = 0:L % degree
%   for m = -l:l % order    
%     if m > 0
%       part = @(x) sqrt(2)*(-1)^m*real(x);
%     elseif m < 0
%       part = @(x) sqrt(2)*(-1)^m*imag(x);
%     else
%       part = @(x) x;
%     end
%     idx = l*(l+1)+1+m;
%     m = abs(m);
% 
%     Pml = legendre(l,cos(EL + pi/2)); % Associated Legendre Polynomials of degree l
%     if l~=0, Pml = Pml(m + 1,:); end % pick the order m term
% 
%     C = sqrt((2*l + 1)/(4*pi) * gamma(l - m + 1)/gamma(l + m + 1)); % normalization constant
%     Yml = C*Pml'.*exp(1i*m*AZ); % complex spherical harmonics
%     realY(:,idx) = part(Yml); % convert to real SH basis
%   end

    Pml = legendre(l,cos(EL + pi/2)); % Associated Legendre Polynomials of degree l order 0:l

    % normalization constants
    C = sqrt((2*l + 1)/(4*pi) * gamma(l - (0:l) + 1)./gamma(l + (0:l) + 1));
    % complex spherical harmonics
    Yml = repmat(C,numv,1).*Pml'.*exp(1i*repmat(0:l,numv,1).*repmat(AZ,1,l+1));
    
    % convert to real SH basis
    mneg = l*(l+1)+1+(-l:-1);
    m0   = l*(l+1)+1;
    mpos = l*(l+1)+1+(1:l);
    realY(:,mneg) = repmat(sqrt(2)*(-1).^(-l:-1),numv,1).*imag(Yml(:,(l+1):-1:2));
    realY(:,m0) = Yml(:,1);
    realY(:,mpos) = repmat(sqrt(2)*(-1).^(1:l),numv,1).*real(Yml(:,2:l+1));
    
    % record degree and order for ease of access later
    if nargout > 1
      LM = [LM; repmat(l,2*l+1,1), (-l:l)'];
    end
end
end