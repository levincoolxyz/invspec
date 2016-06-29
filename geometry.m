classdef geometry < handle
   properties
     v,f
     numv,numf
     M,L
     D,V
     sparselim % above which switch to sparse matrices
   end
   methods
      function obj = geometry(v,f)
        obj.v = v;
        obj.f = f;
        obj.numv = size(obj.v,1);
        obj.numf = size(obj.f,1);
        obj = obj.lapbel(obj);
      end
      function obj = lapbel(obj)
        % computes the Laplace-Beltrami operator of discrete surfaces

        % cot of angle at x and dual area associated to x in each face 
        cotan = @(x,y,z) (y-x)'*(z-x)./norm(cross(y-x,z-x),2);
        % dualA = @(x,y,z) (((y-x)'*(y-x))*cotan(z,x,y) + ((z-x)'*(z-x))*cotan(y,z,x))/8;
        dualA = @(x,y,z) norm(cross(y-x,z-x),2)/6; % barycentric

        % construction (naive, too many calls to @cotan ?)
        if obj.numv > obj.sparselim
          obj.M = spalloc(obj.numv,obj.numv,obj.numv);
        else
          obj.M = zeros(obj.numv);
        end
        if obj.numv > obj.sparselim
          obj.L = spalloc(obj.numv,obj.numv,7*obj.numv);
        else
          obj.L = zeros(obj.numv);
        end
        for fi=1:obj.numf
          for idx = 0:2
            i = obj.f(fi,idx+1);
            j = obj.f(fi,mod(idx+1,3)+1);
            k = obj.f(fi,mod(idx+2,3)+1);
            obj.L(i,k) = obj.L(i,k) + .5*cotan(obj.v(j,:)',obj.v(k,:)',obj.v(i,:)');
            obj.L(i,i) = obj.L(i,i) - .5*cotan(obj.v(j,:)',obj.v(k,:)',obj.v(i,:)');
            obj.L(i,j) = obj.L(i,j) + .5*cotan(obj.v(k,:)',obj.v(i,:)',obj.v(j,:)');
            obj.L(i,i) = obj.L(i,i) - .5*cotan(obj.v(k,:)',obj.v(i,:)',obj.v(j,:)');
            obj.M(i,i) = obj.M(i,i) + dualA(obj.v(i,:)',obj.v(j,:)',obj.v(k,:)');
          end
        end
      end
      function obj = spectra(obj,numeig)
        if nargin < 3 || isempty(numeig), numeig = obj.numv; end
        [obj.V,obj.D] = eigvf(obj.L,obj.M,numeig);
      end
      function [obj,s] = refine(obj,s,sthreshold)
        % refine near vertices of high conformal factor
        for fi = 1:size(obj.f,1)
          i = obj.f(fi,:);
          if max(abs(s(i))) > sthreshold % face contains vertex of high conformal factor 
            bpos = sum(obj.v(i,:),1)/3; % barycenter position of said face
            bpos = bpos/norm(bpos); % project it back to sphere
            bs = sum(s(i))/3; % average/interp conformal factor at barycenter
            vtmp = [obj.v; bpos];
            s = [obj.s; bs];
          end
        end
        obj.f = fliplr(convhulln(vtmp));
        obj.v = vtmp((unique(obj.f(:))),:);
        s = s((unique(obj.f(:))),:);
        obj.f = (convhulln(obj.v));
        obj.numv = size(obj.v,1);
        obj.numf = size(obj.f,1);
        obj = obj.lapbel(obj);
      end
      function varargout = eigencost(s,obj,lambda_T,numeig,reg,eig0)
        % function [J,GJ] = eigencost(s,M,L,lambda_T,numeig,reg,eig0)
        % 
        % INPUT
        % s        - (numv x 1) vector of per-vertex conformal factors
        % M, L     - (numv x numv) mass and cotan matrices
        % lambda_T - (numv x 1) target discrete spectrum
        % numeig   - number of eigenvalues considered
        % reg      - coefficient for the regularization term
        % eig0     - the starting eigenvalue being considered
        % OUTPUT
        % J        - cost
        % GJ       - gradient
        % 

        if nargin<5 || isempty(reg), reg = 0; end
        if nargin<6 || isempty(eig0) || eig0<2, eig0 = 2; end
        nums = numel(s);
        %% get eigenvalues and eigenvectors
        % [Vec,lambda] = eigvf(obj.L,diag(1./s)*obj.M,numeig);
        [Vec,lambda] = eigvf(obj.L,diag(1./exp(s))*obj.M,numeig); % log-conformal factors
        %% spectral cost and gradient
        lambda_T = lambda_T(end-numeig+1:end); % chop unused target values
        lambda_diff = lambda-lambda_T; % take the normal difference
        lambda_diff = lambda_diff./lambda_T; % change to relative difference
        % lambda_diff = 1./lambda-1./lambda_T; % take inverse difference

        J = .5*sum(lambda_diff(1:(end-1)).^2); % compute norm squared cost

        LLs = obj.L*obj.L*s;
        J = J + s'*LLs/2*reg;  % apply bi-laplacian regularization
        if (nargout <= 1)
          varargout{1} = J;
        else
          varargout{1} = J;
          GJ = zeros(nums,1);
          m = diag(obj.M);
          parfor j = 1:nums % iterate relevant conformal factor
            sj = s(j);
            Mjj = m(j);
            for i = 1:(numeig-eig0+1) % iterate relevant eigenvalue
              vij = Vec(j,i);
        %       wij = lambda_diff(i)*lambda(i)*vij^2/sj^2*Mjj;   % normal difference
              wij = lambda_diff(i)*lambda(i)*vij^2*exp(-sj)*Mjj; % normal difference (log)
              wij = wij/lambda_T(i);                             % relative difference
              GJ(j) = GJ(j) + wij;
            end
            GJ(j) = GJ(j) + LLs(j)*reg;
          end
          varargout{2} = GJ;
        end
      end
   end
end