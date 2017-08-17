function [v,v_T,v_end,f,f_T,s_end,s_T,J_hist,Jc_hist,...
  D_0,D_T,D_endp,D_end] = main(init_data,target_data,...
  method,imax,aC,bC,tC,etolC,aS,bS,tS,etolS,...
  numeigI,pert,reg,refctl)
% function [v,v_T,v_end,f,f_T,s_end,s_T,J_hist,Jc_hist,...
%   D_0,D_T,D_endp,D_end] = main(init_data,target_data,...
%   method,imax,aC,bC,tC,etolC,aS,bS,tS,etolS,...
%   numeig,pert,reg,refctl)
%
% for help see comments in test_script.m

vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
%% starting input mesh
if init_data.num ~= 4
  % import wavefront object file
  if init_data.num == 1
    fid = fopen(['../meshes/' init_data.dat '.obj'],'rt');
    [v,f] = readwfobj(fid);

  % sphere of chosen size
  elseif init_data.num == 2
    if exist(['../meshes/' init_data.dat '.mat'],'file')
      load(['../meshes/' init_data.dat '.mat']);
    else
      ssize = str2double(init_data.dat);
      v = ParticleSampleSphere('Vo',RandSampleSphere(ssize));
      f = fliplr(convhulln(v));
      save(['../meshes/' init_data.dat '.mat'],'f','v');
    end

  % load face-vertex from *.mat [need v and f]
  elseif init_data.num == 3
    load(['../meshes/' init_data.dat]);
  end
  
  v = v - repmat(volCenter(v,f),size(v,1),1);
  v = v*makeUnitArea(v,f)*sqrt(4*pi);
  [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);
end
%% target spectrum (and shapes)
if target_data.num == 5
  v_T = []; f_T = []; s_T = [];
  D_T = target_data.D_T;
  if numel(D_T) >= numeigI
    D_T = D_T(1:numeigI);
  else
    D_T = [D_0(1:(numeigI-numel(D_T))); D_T];
  end
else
  % perturb with random conformal factors at vertices
  if target_data.num == 1
    s_T = exp(-rand(numv,1)*pert);
    f_T = f;
    conf_T = sqrt(kron(1./s_T',1./s_T)); % averaing conformal factors at vertices to edges
    elsq_T = elsq0.*conf_T(isedge); % apply to linearly indexed edge lengths
    [~,v_Thist] = gradescent(@conformalcost,imax,aC,bC,tC,etolC,0,...
      reshape(v',[],1),isedge,elsq_T);
    v_T = reshape(v_Thist(:,end),3,[])';
    
    v_T = v_T - repmat(volCenter(v_T,f_T),size(v_T,1),1);
    v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);

  % perturb with given scalar field
  elseif target_data.num == 2
    % compute mean curvature vertex normal
    invM = diag(1./diag(M));
    Hn = .5*[invM*L*v(:,1) invM*L*v(:,2) invM*L*v(:,3)];
    H = vnorm(Hn);
    vn = Hn./repmat(H,1,3);
    v_T = v - repmat(target_data.dat(v),1,3).*vn*pert;
    f_T = f;
    v_T = v_T - repmat(volCenter(v_T,f_T),size(v_T,1),1);
    v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);
    [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');
    [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);

  % import wavefront object file
  elseif target_data.num == 3
    if exist(['mcf/' target_data.dat '.mat'],'file')
      if init_data.num == 4
        load(['mcf/' target_data.dat '.mat'],'v_T','f_T','s_T','v');
        v_T = v_T - repmat(volCenter(v_T,f_T),size(v_T,1),1);
        v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);
        f = f_T;
        [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);
        vmcf = v;
      else
        vtmp = v;
        load(['mcf/' target_data.dat '.mat'],'v_T','f_T','s_T');
        v_T = v_T - repmat(volCenter(v_T,f_T),size(v_T,1),1);
        v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);
        vmcf = v;
        v = vtmp;
      end
    else
      fid = fopen(['../meshes/' target_data.dat '.obj'],'rt');
      [v_T,f_T] = readwfobj(fid);
      v_T = v_T - repmat(volCenter(v_T,f_T),size(v_T,1),1);
      v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);
      if init_data.num == 4
        [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');
        save(['mcf/' target_data.dat '.mat'],'v_T','f_T','s_T','v');
        f = f_T;
        [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);
        vmcf = v;
      else
        [s_T,vmcf] = meancurvflow(v_T,f_T,1e5,'c');
      end
    end

  % load face-vertex from *.mat [need v_T and f_T]
  elseif target_data.num == 4
    load(['../meshes/' target_data.dat]);
    v_T = v_T - repmat(volCenter(v_T,f_T),size(v_T,1),1);
    v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);
    if init_data.num == 4
      [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');
      f = f_T;
      [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);
      vmcf = v;
    else
      [s_T,vmcf] = meancurvflow(v_T,f_T,1e5,'c');
    end
  end
  
  [M_T,L_T] = lapbel(v_T,f_T);
  D_T = eigvf(L_T,M_T,numeig);
%   D_T = eigvf(L,diag(1./s_T)*M,numeig); % cheat with cMCF spectra (if init_data=4)
end
%% refinement criterion
if init_data.num == 4 || isempty(refctl)
  scheck = Inf; % disable refinement routines for cMCF'ed mesh
  sthreshold = 0; % abs(log(1/(conformal factors)))
  refinestop = @(x,optimValues,state) max(exp(x))>scheck;
  refineIter = 0;
else
  scheck = refctl(1); % abs(log(1/(conformal factors)))
  sthreshold = refctl(3);
  sgthreshold = refctl(4);
  refinestop = @(x,optimValues,state) max(exp(x))>scheck;
  refineIter = 0;
  maxRefine = refctl(2);
end
%% MIEP2 via gradient / BFGS descent
% s0 = exp(-zeros(numv,1));
s0 = zeros(numv,1); % if using log conformal factors
vlim = max(max(abs(v)));
vlim = [-vlim vlim];
for neig = [numeig]%(unique(round(logspace(log10(2),log10(numeig),5))))%[2:20 40:20:numeig]
  if neig < .2*numv, reg = 0; end
  if strcmp(method, 'BFGS')
    test = @(s) eigencost(s,M,L,D_T,neig,reg);
    options = optimset('GradObj','on','display','iter-detailed',...
      'maxiter',imax,'tolFun',etolS,'tolx',etolS,'largescale','off',...
      'outputfcn',refinestop);
    [s,J_hist,exitflag] = fminunc(test,s0,options);
    skipnow = 0;
    while exitflag == -1 && refineIter < maxRefine
      oldv = v;
      [v,f,s] = refine(v,f,s,sthreshold,sgthreshold);
      if numv == size(v,1)
        if norm(oldv-v) < 10*eps
          skipnow = 1;
        end
      end
      if skipnow
          s0 = s;
      else
        trisurf(f,v(:,1),v(:,2),v(:,3),s,...
          'facecolor','interp');
        legend(num2str(size(v,1),'#vtx %d'),'location','best');
        skipstr = input('mesh refined due to bad conformal factors, continue? Y/N [Y]: \n','s');
        
        if ~isempty(target_data.adapt) && target_data.num == 3
          [~,sizefit] = min(abs(target_data.adapt - size(v,1)));
          adaptedTarget = [regexprep(target_data.dat,'\d*','') ...
            num2str(target_data.adapt(sizefit),'%d')];
          fid = fopen(['../meshes/' adaptedTarget '.obj'],'rt');
          [v_T,f_T] = readwfobj(fid);
          [M_T,L_T] = lapbel(v_T,f_T);
        end
        
        if ~isempty(skipstr)
          if skipstr == 'N' || skipstr == 'n'
            options = optimset('GradObj','on','display','iter-detailed',...
              'maxiter',imax,'tolFun',etolS,'tolx',etolS,'largescale','off');
            [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);
            D_T = eigvf(L_T,M_T,numeig);
            test = @(s) eigencost(s,M,L,D_T,numeig,reg);
            [s,J_hist] = fminunc(test,s,options);
            break;
          end
        end
        
        [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);
%         s0 = zeros(numv,1); % if using log conformal factors and restarting
        s0 = s; % if not restarting
        D_T = eigvf(L_T,M_T,numeig);
        test = @(s) eigencost(s,M,L,D_T,numeig,reg);
      end
      [s,~,exitflag] = fminunc(test,s0,options);
      refineIter = refineIter + 1;
    end
    if exitflag ~= 1
      options = optimset('GradObj','on','display','iter-detailed',...
        'maxiter',imax,'tolFun',etolS,'tolx',etolS,'largescale','off');
      [s,J_hist] = fminunc(test,s,options);
    end
  elseif strcmp(method, 'GD')
    [J_hist,s] = gradescent(@eigencost,imax,aS,bS,tS,etolS,0,...
      s0,[],M,L,D_T,neig,reg);
  else
    error('unknown descent method');
  end

  s0 = s(:,end);  

  h = figure(); % manual inspection of conformal factors

  crange = [min([s_T;s0]) max([s_T;s0])];
  subplot(1,2,1); hold all; view(3); grid on; axis equal
  trisurf(f_T,vmcf(:,1),vmcf(:,2),vmcf(:,3),s_T,...
    'facecolor','interp','edgecolor','none');
  set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
  xlabel('x'); ylabel('y'); zlabel('z');
  title('spherical/cMCF mesh');
  caxis(crange);
  colorbar('southoutside')

  subplot(1,2,2); hold all; view(3); grid on; axis equal
  trisurf(f,v(:,1),v(:,2),v(:,3),s0,...
    'facecolor','interp','edgecolor','none');
  set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
  xlabel('x'); ylabel('y'); zlabel('z');
  title('spectrally opt. mesh');
  caxis(crange);
  colorbar('southoutside')
  
  colormap jet
  pause(.1);
%   s_err = norm(s_T - s0)./norm(s_T)*100;
%   fprintf('total conformal factor error = %g\n',s_err);
  if isempty(refctl)
    skipstr = 'y';
  else
    skipstr = input('Continue optimization? Y/N [Y]: \n','s');
  end
  close(h);
  if ~isempty(skipstr), if skipstr == 'N' || skipstr == 'n', break; end; end
end
% s_end = s(:,end);
s_end = exp(s(:,end)); % if using log conformal factors
D_endp = eigvf(L,diag(1./s_end)*M,numeig);
%% record conformal factor changes during descent (if using in-house algorithm)
% figure(); view(3); grid on; axis equal
% h = trisurf(f,v(:,1),v(:,2),v(:,3),s(:,1),...
%   'facecolor','interp','edgecolor','none');
% vlim = max(max(abs(v)));
% vlim = [-vlim vlim];
% hh = gca;
% set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
% xlabel('x'); ylabel('y'); zlabel('z');
% colormap jet;
% for i=[1:9 unique(round(logspace(log10(10),log10(numel(s,2)),60)))];
% h.CData = s(:,i);
% hh.View = hh.View + [-5 5];
% hgexport(gcf,num2str(i,'descendence%5d.png'),...
%   hgexport('factorystyle'), 'Format', 'png');
% end
% unix('convert -delay 50 -loop 0 descendence*.png descendence.gif');
% unix('rm -f descendence*.png');
%% conformal embedding/fit
conf = sqrt(kron(1./s_end',1./s_end)); % averaing conformal factors at vertices to edges
elsq_end = elsq0.*conf(isedge); % apply to linearly indexed edge lengths
[Mc,Lc,tc] = el2ew(f,elsq_end); % test if averaging conformal factors is good enough
figure(); view(3); axis equal;
trisurf(f,v(:,1),v(:,2),v(:,3),tc);
D_endpp = eigvf(Lc,Mc,numeig);
figure();
hold all;
plot(imag(D_endpp))
plot(-real(D_endpp))
plot(-flipud(real(D_endp)))
legend('imaginary averaged spectrum','real averaged spectrum','pre-averaging spectrum');
if strcmp(method, 'BFGS')
  test = @(v) conformalcost(v,isedge,elsq_end);
  options = optimset('GradObj','on','display','iter-detailed',...
    'maxiter',imax,'tolFun',etolC,'tolx',etolC,'largescale','off');
  [vhist,Jc_hist] = fminunc(test,reshape(v',[],1),options);
elseif strcmp(method, 'GD')
  [Jc_hist,vhist] = gradescent(@conformalcost,imax,aC,bC,tC,etolC,0,...
    reshape(v',[],1),[],isedge,elsq_end);
end
v_end = reshape(vhist(:,end),3,[])';
[M_end,L_end] = lapbel(v_end,f);
D_end = eigvf(L_end,M_end,numeig);

end

% initialize crucial data about a given mesh (really should be made object-oriented)
function [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeig)
%% learn to count
numv = size(v,1); % number of vertices
numf = size(f,1); % number of faces
numeig = numv*(numeig <= 0) + ...
  ceil(numv*numeig)*(numeig <= 1 && numeig > 0) + ...
  min(numv,numeig)*(numeig > 1);
%% when is there an edge (mild redundancy)
isedge = zeros(numv);
for fi = 1:numf
  for idx = 0:2
    i = f(fi,idx+1);
    j = f(fi,mod(idx+1,3)+1);
    isedge(i,j) = 1;
  end
end
isedge = tril(isedge); % reduce redundancy
isedge = find(isedge); % linear indices
%% compute initial edge lengths squared
elsq0 = zeros(numv);
temp = mod(isedge-1,numv)+1;
temp2 = fix((isedge-1)/numv);
for j = 1:numv
  for i = temp(temp2 == (j-1))'
    elsq0(i,j) = sum((v(i,:)-v(j,:)).^2);
  end
end
elsq0 = nonzeros(elsq0);
%% compute laplacian
[M,L] = lapbel(v,f);
try
  D_0 = eigvf(L,M,numeig);
catch
  disp shit eigs
end
end

% refinement subroutine (currently works only for spherical meshes due to use of convhull)
function [v,f,s] = refine(v,f,s,scrit,sgcrit)
% oldfaces = [];
highval = 1;
highgrad = [0 0 0];
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
for fi = 1:size(f,1)
  i = f(fi,1);
  j = f(fi,2);
  k = f(fi,3);
  if max(abs(alldiff(s(i)))) > scrit % face contains vertices with high conformal factors
%   if max(abs((s(i)))) > scrit % face contains vertices with high conformal factors
    bpos = sum(v([i j k],:),1)/3; % barycenter position of said face
%     bpos = bpos/norm(bpos)*mean(vnorm(v([i j k],:))); % project it back to sphere
    bpos = bpos/norm(bpos)*mean(vnorm(v)); % project it back to sphere
    bs = sum(s([i j k]))/3; % average/interp conformal factor at barycenter
    v = [v; bpos];
    s = [s; bs];
    highval = 1;
  end
  if abs(s(i)-s(j)) > sgcrit % edge ij has high conformal factor gradient
    bpos = sum(v([i j],:),1)/2; % midpoint of said edge
%     bpos = bpos/norm(bpos)*mean(vnorm(v([i j],:))); % project it back to sphere
    bpos = bpos/norm(bpos)*mean(vnorm(v)); % project it back to sphere
    bs = sum(s([i j]))/2; % average/interp conformal factor
    v = [v; bpos];
    s = [s; bs];
    highgrad(1) = 1;
  end
  if abs(s(j)-s(k)) > sgcrit % edge jk has high conformal factor gradient
    bpos = sum(v([j k],:),1)/2; % midpoint of said edge
%     bpos = bpos/norm(bpos)*mean(vnorm(v([j k],:))); % project it back to sphere
    bpos = bpos/norm(bpos)*mean(vnorm(v)); % project it back to sphere
    bs = sum(s([j k]))/2; % average/interp conformal factor
    v = [v; bpos];
    s = [s; bs];
    highgrad(2) = 1;
  end
  if abs(s(k)-s(i)) > sgcrit % edge ki has high conformal factor gradient
    bpos = sum(v([k i],:),1)/2; % midpoint of said edge
%     bpos = bpos/norm(bpos)*mean(vnorm(v([k i],:))); % project it back to sphere
    bpos = bpos/norm(bpos)*mean(vnorm(v)); % project it back to sphere
    bs = sum(s([k i]))/2; % average/interp conformal factor
    v = [v; bpos];
    s = [s; bs];
    highgrad(3) = 1;
  end
  
%   if highval && sum(highgrad) == 3
%     oldfaces = [oldfaces; fi];
%     nv = size(v,1);
%     f = [f; i nv-2 nv-3;
%       nv-2 j nv-3;
%       j nv-1 nv-3;
%       nv-1 k nv-3;
%       k nv nv-3;
%       nv i nv-3];
%   elseif ~highval && sum(highgrad) == 3
%     oldfaces = [oldfaces; fi];
%     nv = size(v,1);
%     f = [f; i nv-2 nv;
%       nv-2 j nv-1;
%       nv-1 k nv;
%       nv nv-2 nv-1];
%   end
end
% f(oldfaces,:) = [];

% in case mesh is 'spherical'
f = fliplr(convhulln(v));
v = v(sort(unique(f(:))),:);
s = s(sort(unique(f(:))),:);
f = fliplr(convhulln(v));

% DT = delaunayTriangulation(v);
% f = DT.ConnectivityList;
end

% refinement wrapper for in-house gradient descent !INCOMPLETE!
function [stop,s,M,L,D_T,neig,reg] = refineGD(s,M,L,D_T,neig,reg)
  stop = 0;
end

% calculates all possible differences of elements in a vector
function [xdiff] = alldiff(x)
x = x';
xdiff = arrayfun(@(k) x(k:end)-x(k), 1:numel(x),'un',0);
xdiff = [xdiff{:}];
xdiff = xdiff';
end