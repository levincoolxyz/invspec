function [v,v_T,v_end,f,f_T,a_end,s_end,s_T,J_hist,Jc_hist,...
  D_0,D_T,D_endp,D_end,vmcf] = mainSH(init_data,target_data,...
  method,imax,aC,bC,tC,etolC,aS,bS,tS,etolS,...
  numeigI,maxL,numa,pert)
% function [v,v_T,v_end,f,f_T,s_end,s_T,J_hist,Jc_hist,...
%   D_0,D_T,D_endp,D_end] = mainSH(init_data,target_data,...
%   method,imax,aC,bC,tC,etolC,aS,bS,tS,etolS,...
%   numeigI,pert,maxL)
%
% for help see comments in test_scriptSH.m

vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
numL = (maxL+1)^2;
mcf_step = 1e5;
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
  Y_v = sphericalHarmonicBase(v,maxL);
end
%% target spectrum (and shapes)
if target_data.num == 5
  v_T = []; f_T = []; a_T = target_data.a;
  D_T = target_data.dat;
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
    [s_T,vmcf] = meancurvflow(v_T,f_T,mcf_step,'c');
%     [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);

  % import wavefront object file
  elseif target_data.num == 3
    if exist(['mcf/' target_data.dat '.mat'],'file')
      if init_data.num == 4
        load(['mcf/' target_data.dat '.mat'],'v_T','f_T','s_T','v','sa_T');
        v_T = v_T - repmat(volCenter(v_T,f_T),size(v_T,1),1);
        v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);
        f = f_T;
        vmcf = v;
        [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);
        Y_v = sphericalHarmonicBase(v,maxL);
      else
        vtmp = v;
        load(['mcf/' target_data.dat '.mat'],'v_T','f_T','s_T','v','sa_T');
        v_T = v_T - repmat(volCenter(v_T,f_T),size(v_T,1),1);
        v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);
        vmcf = v;
        v = vtmp;
      end
    else
      fid = fopen(['../meshes/' target_data.dat '.obj'],'rt');
      [v_T,f_T] = readwfobj(fid);
      v_T = v_T - repmat(volCenter(v_T,f_T),numv,1);
      v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);
      if init_data.num == 4
        [s_T,v] = meancurvflow(v_T,f_T,mcf_step,'c');
        save(['mcf/' target_data.dat '.mat'],'v_T','f_T','s_T','v','sa_T');
        f = f_T;
        vmcf = v;
        [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);
        Y_v = sphericalHarmonicBase(v,maxL);
      else
        [s_T,vmcf] = meancurvflow(v_T,f_T,mcf_step,'c');
      end
    end

  % load face-vertex from *.mat [need v_T and f_T]
  elseif target_data.num == 4
    load(['../meshes/' target_data.dat]);
    v_T = v_T - repmat(volCenter(v_T,f_T),size(v_T,1),1);
    v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);
    if init_data.num == 4
      [s_T,v] = meancurvflow(v_T,f_T,mcf_step,'c');
      f = f_T;
      vmcf = v;
      [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeigI);
      Y_v = sphericalHarmonicBase(v,maxL);
    else
      [s_T,vmcf] = meancurvflow(v_T,f_T,mcf_step,'c');
    end
  end
  
%   M_T = lapbel(vmcf,f_T);
%   Y_vmcf = sphericalHarmonicBase(vmcf,maxL);
%   delta = Y_vmcf'*M_T*Y_vmcf;
%   a_T = delta\(Y_vmcf'*M_T*s_T);
%   D_T = eigvfSH(a_T,numeig,maxL);
  
  [M_T,L_T] = lapbel(v_T,f_T);
  D_T = eigvf(L_T,M_T,numeig);
end
% s_T = Y_v*a_T;
%% MIEP2 via gradient / BFGS descent
% a0 = [2*sqrt(pi);zeros(numL-1,1)];
% a0 = [2*sqrt(pi);ones(numL-1,1)./(2:numL)'];
% a0 = [2*sqrt(pi);ones(numL-1,1)./(2:numL).^2'];
% a0 = [2*sqrt(pi);1e-3*ones(numL-1,1)./(2:numL)'];
% a0 = [2*sqrt(pi);ones(numeig-1,1)./exp(2:numeig)';zeros(numL-numeig,1)];
% a0 = [2*sqrt(pi);1e-3*ones(numeig-1,1)./exp(2:numeig)';zeros(numL-numeig,1)];
% a0 = [2*sqrt(pi);1e-3*ones(numeig-1,1)./exp(2:numeig).^2';zeros(numL-numeig,1)];
% a0 = [2*sqrt(pi);eps*ones(numL-1,1)];
a0 = [2*sqrt(pi);eps*ones(numa-1,1)];

if strcmp(method, 'BFGS')
  test = @(a) eigencostSH(a,D_T,numeig,maxL);
  options = optimset('GradObj','on','display','iter-detailed',...
    'maxiter',imax,'tolFun',etolS,'tolx',etolS,'largescale','off');
  [a,J_hist] = fminunc(test,a0,options);
elseif strcmp(method, 'GD')
  [J_hist,a] = gradescent(@eigencostSH,imax,aS,bS,tS,etolS,0,...
    a0,[],D_T,numeig,maxL);
else
  error('unknown descent method');
end

a_end = [a(:,end);zeros(numL-numa,1)];
D_endp = eigvfSH(a_end,numeig,maxL);

s_end = Y_v*a_end;
%% conformal embedding/fit
conf = sqrt(kron(1./s_end',1./s_end)); % averaing conformal factors at vertices to edges
elsq_end = elsq0.*conf(isedge); % apply to linearly indexed edge lengths
[Mc,Lc,tc] = el2ew(f,elsq_end); % test if averaging conformal factors is good enough
D_endpp = eigvf(Lc,Mc,numeig);
if norm(imag(D_endpp)) >= 10*eps % pathological averaging
  warning('pathological averaging (nontrivial imaginary spectrum)');
  figure();
  subplot(1,2,1); hold all; view(3); axis equal;
  trisurf(f,v(:,1),v(:,2),v(:,3),tc);
  subplot(1,2,2); hold all; grid on;
  plot(imag(D_endpp))
  plot(-real(D_endpp))
  plot(-flipud(real(D_endp)))
  legend('imaginary averaged spectrum','real averaged spectrum','pre-averaging spectrum',...
    'location','best');
  D_end = [];
  v_end = [];
  Jc_hist = Inf;
  pause(1);
else
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

v_T = v_T - repmat(volCenter(v_T,f_T),size(v_T,1),1);
v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);
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