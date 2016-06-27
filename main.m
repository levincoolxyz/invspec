function [v,v_T,v_end,f,f_T,s_end,s_T,J_hist,Jc_hist,...
  D_0,D_T,D_endp,D_end] = main(init_data,target_data,...
  imax,aC,bC,tC,etolC,aS,bS,tS,etolS,...
  numeig,pert,reg)
% function [v,v_T,v_end,f,f_T,s_end,s_T,J_hist,Jc_hist,...
%   D_0,D_T,D_endp,D_end] = main(init_data,target_data,...
%   imax,aC,bC,tC,etolC,aS,bS,tS,etolS,...
%   numeig,pert,reg)

vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
%% initial mesh
if target_data.num ~= 3 && target_data.num ~= 4
  % import wavefront object file
  if init_data.num == 1
    fid = fopen(['../meshes/' init_data.dat '.obj'],'rt');
    [v,f] = readwfobj(fid);

  % sphere of chosen size
  elseif init_data.num == 2
    if exist(['../meshes/' init_data.dat '.mat'],'file')
      load(['../meshes/' init_data.dat '.mat']);
    else
      ssize = str2num(init_data.dat);
      v = ParticleSampleSphere('Vo',RandSampleSphere(ssize));
      f = fliplr(convhulln(v));
      save(['../meshes/' init_data.dat '.mat'],'f','v');
    end

  % load face-vertex from *.mat
  elseif init_data.num == 3
    load(['../meshes/' init_data.dat]);
  end
  [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeig);
end
%% target spectrum (+mesh for testing)
if isfield(target_data,'D_T')
  v_T = []; f_T = []; s_T = [];
  D_T = target_data.D_T;
  if numel(D_T) >= numeig
    D_T = D_T(1:numeig);
  else
    D_T = [D_0(1:(numeig-numel(D_T))); D_T];
  end
else
  % perturb with random conformal factors at vertices
  if target_data.num == 1
    s_T = exp(-rand(numv,1)*pert);
    f_T = f;
    conf_T = sqrt(kron(1./s_T',1./s_T));
  %   elsq_T = elsq0.*conf_T;
    elsq_T = elsq0.*conf_T(isedge); % linear indices
    [~,v_Thist] = gradescent(@conformalcost,imax,aC,bC,tC,etolC,0,...
      reshape(v',[],1),isedge,elsq_T);
    v_T = reshape(v_Thist(:,end),3,[])';

  % perturb with given scalar field
  elseif target_data.num == 2
    % compute mean curvature vertex normal
    Hn = .5*[inv(M)*L*v(:,1) inv(M)*L*v(:,2) inv(M)*L*v(:,3)];
    H = vnorm(Hn);
    vn = Hn./repmat(H,1,3);
    v_T = v - repmat(target_data.dat(v),1,3).*vn*pert;
    f_T = f;
    [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');
    [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeig);

  % import wavefront object file
  elseif target_data.num == 3
    if exist(['mcf/' target_data.dat '.mat'],'file')
      load(['mcf/' target_data.dat '.mat']);
    else
      fid = fopen(['../meshes/' target_data.dat '.obj'],'rt');
      [v_T,f_T] = readwfobj(fid);
      [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');
      save(['mcf/' target_data.dat '.mat'],'v_T','f_T','s_T','v');
    end
    f = f_T;
    [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeig);

  % load face-vertex from *.mat
  elseif target_data.num == 4
    load(['../meshes/' target_data.dat]);
    [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');
    f = f_T;
    [numv,numeig,isedge,elsq0,M,L,D_0] = initialize(v,f,numeig);
  end
  
  [M_T,L_T] = lapbel(v_T,f_T);
  D_T = eigvf(L_T,M_T,numeig);

%   D_T = eigvf(L,diag(1./s_T)*M,numeig);
end
%% MIEP2 via naive gradient descent
% s0 = exp(-zeros(numv,1));
s0 = zeros(numv,1); % if using log-conformal factors
vlim = max(max(abs(v)));
vlim = [-vlim vlim];
for neig = [numeig 2:2]%(unique(round(logspace(log10(2),log10(numeig),5))))%[2:20 40:20:numeig]
  if neig < 20, reg = 0; end
  test = @(s) eigencost(s,M,L,D_T,neig,reg);
  options = optimset('GradObj','on','display','iter-detailed',...
    'maxiter',imax,'tolFun',etolS,'tolx',etolS,'largescale','off');
  [s,J_hist] = fminunc(test,s0,options);

  % [J_hist,s] = gradescent(@eigencost,imax,aS,bS,tS,etolS,0,...
  %   s0,M,L,D_T,neig,reg);

  % s = s_T; J_hist = [];
  s0 = s(:,end);  

  figure(1); % manual inspection of conformal factors

  crange = [min([s_T;s0]) max([s_T;s0])];
  subplot(1,2,1); hold all; view(3); grid on; axis equal
  trisurf(f,v(:,1),v(:,2),v(:,3),s_T,...
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
  skipstr = input('Continue spectral optimization? Y/N [Y]: \n','s');
  close(1);
  if ~isempty(skipstr), if skipstr == 'N' || skipstr == 'n', break; end; end
end
% s_end = s(:,end);
s_end = exp(s(:,end)); % if using log-conformal factors
D_endp = eigvf(L,diag(1./s_end)*M,numeig);
%% record descent conformal changes
% figure(); view(3); grid on; axis equal
% h = trisurf(f,v(:,1),v(:,2),v(:,3),s(:,1),...
%   'facecolor','interp','edgecolor','none');
% vlim = max(max(abs(v)));
% vlim = [-vlim vlim];
% hh = gca;
% set(gca,'xlim',vlim,'ylim',vlim,'zlim',vlim);
% xlabel('x'); ylabel('y'); zlabel('z');
% colormap jet;
% for i=[1:9 10:10:100 200:100:numel(s,2)];
% h.CData = s(:,i);
% hh.View = hh.View + [-5 5];
% hgexport(gcf,num2str(i,'descendence%5d.png'),...
%   hgexport('factorystyle'), 'Format', 'png');
% end
% unix('convert -delay 50 -loop 0 descendence*.png descendence.gif');
% unix('rm -f descendence*.png');
%% debug error checking
% figure();plot(s_T);hold all;plot(s0);plot(s_end,'x')
% figure();plot(D_T);hold all;plot(D_0);plot(D_endp,'x')
% s_err = norm(s_T - s_end)./norm(s_T)*100;
% fprintf('conformal factors error = %g%%\n',s_err);
%% conformal embedding/fit
conf = sqrt(kron(1./s_end',1./s_end));
% elsq_end = elsq0.*conf;
elsq_end = elsq0.*conf(isedge); % linear indices
% elsq_end = elsq0; % when testing / don't know if above converges in a nice fashion
test = @(v) conformalcost(v,isedge,elsq_end);
options = optimset('GradObj','on','display','iter-detailed',...
  'maxiter',imax,'tolFun',etolC,'tolx',etolC,'largescale','off');
[vhist,Jc_hist] = fminunc(test,reshape(v',[],1),options);
% [Jc_hist,vhist] = gradescent(@conformalcost,imax,aC,bC,tC,etolC,0,...
%   reshape(v',[],1),isedge,elsq_end);
v_end = reshape(vhist(:,end),3,[])';
[M_end,L_end] = lapbel(v_end,f);
D_end = eigvf(L_end,M_end,numeig);

end

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
D_0 = eigvf(L,M,numeig);
end