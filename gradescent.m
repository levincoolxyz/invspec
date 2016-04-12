function [J,v] = gradescent(costf,imax,alpha,beta,v0,varargin)
% [J,v] = gradescent(costf,imax,alpha,beta,v0,varargin)
% gradient descent with backtracking
% 
% INPUTS
% costf       - string or handle of the cost function which returns a
%               scalar cost and a gradient (nx1) vector in that order
% imax        - the maximum allowed number of iteration
% alpha,beta  - backtracking control parameters (cf. Section 5.1.2 in
%               https://www.cs.cmu.edu/~ggordon/10725-F12/scribes/10725_Lecture5.pdf)
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
for i = 1:imax
  [J(i),GJ] = feval(costf,v(:,i),varargin{:});
  if i>10
    if J(i) - J(i-1) == 0
      J(i:end) = [];
      v(:,i:end) = [];
      break;
    end
  end
  t = (prev_ln>10)*beta^(prev_ln-3) + (prev_ln<=10);
%   t = 1;
  vnext = v(:,i) - t*GJ;
  while feval(costf,vnext,varargin{:}) > J(i) - alpha*t*sum(GJ.^2)
    t = beta*t;
    vnext = v(:,i) - t*GJ;
  end
  v(:,i+1) = vnext;
  prev_ln = ceil(log(t)/log(beta))+1; %previous ln search, use as heuristic
  fprintf('descent iter#%d; J = %g; step = beta^%d\n',i,J(i),prev_ln);
end
