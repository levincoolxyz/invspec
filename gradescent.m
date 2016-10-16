function [J,v] = gradescent(costf,imax,c1,c2,t0,...
  etol,figfg,v0,outputFcn,varargin)
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
% outputFcn   - output function called at the end of each step
% varargin    - additional parameters passed on to costf
% 
% OUTPUTS
% J           - (imax x 1) vector of the optimized cost history
% v           - (n x imax) matrix of the optimizing parameter history
%

% c3 = .5;
% c4 = 2;
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
  	if abs((J(i-1) - J(i))/J(i-1)) <= etol || norm(GJ) <=eps
      fprintf('Converged at iter#%d; J = %g; |GJ| = %g; tau = %g\n',...
        i,J(i),norm(GJ),tau);
      i = i + (i < imax);
      J(i:end) = [];
      v(:,i:end) = [];
      break;
    end
  end
  while abs((tau-tau_old)/tau)>eps && abs(tau) > eps
    if abs(tau) < 1e-10
      warning('wee:too:low','small stepsize. tau = %g',tau);
%       tau = 1;
      break;
%    elseif abs((tau - tau_old)/tau) > 1e-5
%       fprintf('try tau = %g\n',tau);
%    else
%      fprintf('.');
    end
    tau_old = tau;
    vp1 = v(:,i) + tau*u;
    [Jp1,GJp1] = feval(costf,vp1,varargin{:});
    if Jp1 > J(i) + c1*tau*dot(GJ,u)
      b = tau;
    elseif dot(GJp1,u) < c2*dot(GJ,u)
      a = tau;
    else
      break;
    end
    if b < inf
%       tau = (1-c3)*a+c3*b;
      tau = (a+b)/2;
    else
%       tau = a*c4;
      tau = 2*a;
    end
  end
%   tau = 1e-3;
  fprintf('\n');
  v(:,i+1) = v(:,i) + tau*u;
  fprintf('descent iter#%d; J = %g; |GJ| = %g; tau = %g\n',...
    i,J(i),norm(GJ),tau);
  if ~isempty(outputFcn)
    [stop,vnew,varargin] = feval(outputFcn,v(:,i+1),varargin);
    if numel(vnew) > size(v,1)
      v = [v;zeros(numel(vnew)-size(v,1),imax)];
      v(:,i+1) = vnew;
    elseif numel(vnew) < size(v,1)
      v = v(numel(vnew),:);
      v(:,i+1) = vnew;
    end
    if stop, break; end
  end
end

if(figfg)
  figure(); hold all; grid on;
  plot(J,'k.:')
  set(gca,'xscale','log','yscale','log'); axis square
  % linearregress(log(1:imax),log(J),1);
end
