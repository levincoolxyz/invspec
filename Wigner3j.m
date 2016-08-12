function W = Wigner3j( J, M )
% function W = Wigner3j( J, M ) 
% Compute the Wigner 3j symbol using the Racah formula. 
%
% INPUT
% J - [J1, J2, J3] \in \Z^3
% M - [M1, M2, M3] \in \Z^3[1/2]
% 
% OUTPUT
% W - Wigner 3j symbol
%
% According to seletion rules, W = 0 unless:
%   |Ji - Jj| <= Jk <= (Ji + Jj)    ({i,j,k} are permutations of {1,2,3})
%   |Mi| <= Ji    (i = 1,2,3)
%   M1 + M2 + M3 = 0
% 
% Reference: 
% Wigner 3j-Symbol entry of Eric Weinstein's Mathworld:
% http://mathworld.wolfram.com/Wigner3j-Symbol.html
%
% Inspired by Wigner3j.m by David Terr, Raytheon, 6-17-04
%
% Kobi Kraus, Technion, Updated 1-8-13.

j1 = J(1); j2 = J(2); j3 = J(3);
m1 = M(1); m2 = M(2); m3 = M(3);

% Input error checking
if any( J < 0 ),
    error( 'The j must be non-negative' )
elseif any( rem( [J, M], 0.5 ) ),
    error( 'All arguments must be integers or half-integers' )
elseif any( rem( (J - M), 1 ) )
    error( 'j123 and m123 do not match' );
end

% Selection rules
if ( j3 > (j1 + j2) ) || ( j3 < abs(j1 - j2) ) ... % j3 out of interval
   || ( m1 + m2 + m3 ~= 0 ) ... % non-conserving angular momentum
   || any( abs( M ) > J ), % m is larger than j
    W = 0;
    return
end
    
% Simple common case
if ~any( M ) && rem( sum( J ), 2 ), % m1 = m2 = m3 = 0 & j1 + j2 + j3 is odd
    W = 0;
    return
end

% Evaluation
t1 = j2 - m1 - j3;
t2 = j1 + m2 - j3;
t3 = j1 + j2 - j3;
t4 = j1 - m1;
t5 = j2 + m2;

tmin = max( 0,  max( t1, t2 ) );
tmax = min( t3, min( t4, t5 ) );

t = tmin : tmax;

% W = sum( (-1).^t .* ...
%   exp( - ones(1,6) * gammaln( ([t; t-t1; t-t2; t3-t; t4-t; t5-t]) + 1 ) + ...
%   0.5 * [-1, ones(1,9)] * gammaln(([j1+j2+j3+1; j1+j2-j3; j1-j2+j3; -j1+j2+j3;
%   j1+m1; j1-m1; j2+m2; j2-m2; j3+m3; j3-m3]) + 1) ) ) * (-1)^( j1-j2-m3 );

% digits(20); % vpa precision default 32
W = sum( (-1).^t .* ...
  exp( - ones(1,6) * gammaln( vpa([t; t-t1; t-t2; t3-t; t4-t; t5-t]) + 1 ) + ...
  0.5 * [-1, ones(1,9)] * gammaln(vpa([j1+j2+j3+1; j1+j2-j3; j1-j2+j3; -j1+j2+j3;
  j1+m1; j1-m1; j2+m2; j2-m2; j3+m3; j3-m3]) + 1) ) ) * (-1)^( j1-j2-m3 );

% factl = @(x) gammaln(x+1);
% W = sum( (-1).^t .* ...
%   exp( - factl(t) - factl(t-t1) - factl(t-t2) - factl(t3-t) - factl(t4-t) - factl(t5-t) + ...
%        (- factl(j1+j2+j3+1) + factl(j1+j2-j3) + factl(j1-j2+j3) + factl(-j1+j2+j3) + ...
%        factl(j1+m1) + factl(j1-m1) + factl(j2+m2) + factl(j2-m2) + factl(j3+m3) + factl(j3-m3))...
%        /2)) * (-1)^( j1-j2-m3 );

% fact = @(x) gamma(x+1);
% W = sum( (-1).^t ...
%   ./ fact(t) ./ fact(t-t1) ./ fact(t-t2) ./ fact(t3-t) ./ fact(t4-t) ./ fact(t5-t) ...
%        ./ sqrt(fact(j1+j2+j3+1)) .* sqrt(fact(j1+j2-j3)) .* sqrt(fact(j1-j2+j3)) .* ...
%        sqrt(fact(-j1+j2+j3)) .* sqrt(fact(j1+m1)) .* sqrt(fact(j1-m1)) .* sqrt(fact(j2+m2)) .* ...
%        sqrt(fact(j2-m2)) .* sqrt(fact(j3+m3)) .* sqrt(fact(j3-m3))) * (-1)^( j1-j2-m3 );

% W = sum( (-1).^t .* ...
%   prod(-factorial( vpi([t; t-t1; t-t2; t3-t; t4-t; t5-t]) )) * ...
%        prod(factorial( vpi([j1+j2+j3-1, j1+j2-j3, j1-j2+j3, -j1+j2+j3, ...
%                 j1+m1, j1-m1, j2+m2, j2-m2, j3+m3, j3-m3]) ) ...
%        .*[-1,ones(1,9)] ./ 2 ) ) * (-1)^( j1-j2-m3 );
         
% Warnings
if isnan( W )
    warning( 'MATLAB:Wigner3j:NaN', 'Wigner3J is NaN!' )
elseif isinf( W )
    warning( 'MATLAB:Wigner3j:Inf', 'Wigner3J is Inf!' )
end