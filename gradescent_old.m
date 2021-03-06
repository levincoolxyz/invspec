function [J,v] = gradescent_old(costf,imax,alpha,beta,t0,...
  etol,figfg,v0,varargin)
% [J,v] = gradescent(costf,imax,alpha,beta,t0,...
%   etol,figfg,v0,varargin)
% gradient descent with backtracking
% 
% INPUTS
% costf       - string or handle of the cost function which returns a
%               scalar cost and a gradient (nx1) vector in that order
% imax        - the maximum allowed number of iteration
% alpha,beta  - backtracking control parameters (cf. Section 5.1.2 in
%               https://www.cs.cmu.edu/~ggordon/10725-F12/scribes/10725_Lecture5.pdf)
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

prev_ln = 1;
J = zeros(imax,1);
v = zeros(numel(v0),imax);
v(:,1) = v0;
fprintf('\nAbout to descent, first step might be slow...\n');
for i = 1:imax
  [J(i),GJ] = feval(costf,v(:,i),varargin{:});
  if i>1
  	if (J(i-1) - J(i))/J(i-1) <= etol || (J(i-1) - J(i)) <=10*eps
      fprintf('Converged at iter#%d; J = %g\n\n',i,J(i));
      i = i + (i < imax);
      J(i:end) = [];
      v(:,i:end) = [];
      break;
    end
  end
  t = ((prev_ln>4)*beta^(prev_ln-3) + (prev_ln<=4))*t0;
%   t = t0;
  vnext = v(:,i) - t*GJ;
  prev_ln = 1;
  while feval(costf,vnext,varargin{:}) > J(i) - alpha*t*sum(GJ.^2)
    t = beta*t;
    vnext = v(:,i) - t*GJ;
    prev_ln = prev_ln + 1;
  end
  v(:,i+1) = vnext;
  fprintf('descent iter#%d; J = %g; stepsize *= %g\n',...
    i,J(i),t0*beta^(prev_ln-1));
end

if(figfg)
  figure(); hold all; grid on;
  plot(J,'k.:')
  set(gca,'xscale','log','yscale','log'); axis square
  % linearregress(log(1:imax),log(J),1);
end