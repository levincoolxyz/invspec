function [J,v] = gradescent(costf,imax,c1,c2,t0,...
  etol,figfg,v0,varargin)
% [J,v] = gradescent(costf,imax,c1,c2,t0,...
%   etol,figfg,v0,varargin)
% gradient descent with backtracking
% 
% INPUTS
% costf       - string or handle of the cost function which returns a
%               scalar cost and a gradient (nx1) vector in that order
% imax        - the maximum allowed number of iteration
% c1,c2       - line search control parameters
% t0          - initial stepsize multiplier
% etol        - relative error tolerance
% figfg       - 1 => gives log-log plot of energy, 0 not
% v0          - initial condition (nx1) vector
% varargin    - additional parameters passed on to costf
% 
% OUTPUTS
% J           - (imax x 1) vector of the optimized cost history
% v           - (n x imax) matrix of the optimizing parameter history
%

tau = t0;
J = zeros(imax,1);
v = zeros(numel(v0),imax);
v(:,1) = v0;
fprintf('\nAbout to descent, first step might be slow...\n');
for i = 1:imax
  a = 0;
  b = inf;
  tau_old = 0;
  [J(i),GJ] = feval(costf,v(:,i),varargin{:});
  u = -GJ;
  if i>1
  	if abs((J(i-1) - J(i))/J(i-1)) <= etol || abs(J(i-1) - J(i)) <=10*eps
      fprintf('Converged at iter#%d; J = %g; |GJ| = %g; tau = %g\n',...
        i,J(i),norm(GJ),tau);
      i = i + (i < imax);
      J(i:end) = [];
      v(:,i:end) = [];
      break;
    end
  end
  scale = 2;
  while abs((tau-tau_old)/tau)>eps
    if tau < eps
      fprintf('WARNING: stepsize is getting loooow, tau = %g\n',tau);
    end
    tau_old = tau;
    vp1 = v(:,i) + tau*u;
    [Jp1,GJp1] = feval(costf,vp1,varargin{:});
    up1 = -GJp1;
    if Jp1 > J(i) + c1*tau*dot(GJ,u)
      b = tau;
    elseif dot(GJp1,up1) < c2*dot(GJ,u)
      a = tau;
    else
      break;
    end
    if b < inf
      tau = (a+b)/scale;
    else
      tau = a*scale;
    end
  end
  v(:,i+1) = v(:,i) + tau*u;
  fprintf('descent iter#%d; J = %g; |GJ| = %g; tau = %g\n',...
    i,J(i),norm(GJ),tau);
end

if(figfg)
  figure(); hold all; grid on;
  plot(J,'k.:')
  set(gca,'xscale','log','yscale','log'); axis square
  % linearregress(log(1:imax),log(J),1);
end