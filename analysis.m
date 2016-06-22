clear;
vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
%% control parameters
imax = 5e4; % gradient descent maximum iterations
aC = .7; bC = .8; tC = 10; etolC = 1e-5; % Conformal descent control
aS = .8; bS = .9; tS = 150; etolS = 1e-8; % invSpec descent control
% numeig = .6; % number of eigenvalues used, <=1 => percent, <=0 => all
% pert = .512; % scaling coefficient used to control target perturbation
rng(1432543); % rand seed

%% test perturbation effect on eigenvalues
% close all;
% figure(); hold all; grid on;
% numeig = ceil(.6*200);
% for pert = [0.4 .6 .7 .8 .9 1 1.5]
%   load(num2str(pert,'i2_200_t2_abs(Y32(v))_e0.6p%g.mat'));
%   [M,L] = lapbel(v,f);
%   [M_T,L_T] = lapbel(v_T,f_T);
%   [M_end,L_end] = lapbel(v_end,f);
%   D_T = eigvf(L_T,M_T,numeig);
% 
% %   figure(); hold all; grid on;
% %   ylim([0 .2]);
% %   plot(diag(M))
% %   plot(diag(M_T))
% %   plot(diag(M_end))
%   plot(D_T);
%   pause(1);
% end
% 
%% test mesh refinement laplacian eig convergence
% % close all;
% % figure(); hold all; grid on;
% % vnorm = @(v) sqrt(v(:,3).^2+v(:,1).^2+v(:,2).^2);
% % Y32 = @(v) (v(:,1).^2-v(:,2).^2).*v(:,3)./(vnorm(v)).^3;
% % target_data.dat = @(v) abs(Y32(v));
% % numeig = 100;
% % pert = .6; % scaling coefficient used to control target perturbation
% % for ssize = 1e2:1e2:1e3
% %   v = ParticleSampleSphere('Vo',RandSampleSphere(ssize));
% %   f = fliplr(convhulln(v));
% %   [M,L] = lapbel(v,f);
% %   D(:,ssize/100) = eigvf(L,M,numeig);
% %   Hn = .5*[inv(M)*L*v(:,1) inv(M)*L*v(:,2) inv(M)*L*v(:,3)];
% %   H = vnorm(Hn);
% %   vn = Hn./repmat(H,1,3);
% %   v_T = v - repmat(target_data.dat(v),1,3).*vn*pert;
% %   f_T = f;
% %   [M_T,L_T] = lapbel(v_T,f_T);
% %   D_T = eigvf(L_T,M_T,numeig);
% %   plot(D_T);
% % %   pause(1);
% % end
% % 
% % legend('100','200','300','400','500','600','700','800','900','1E3',...
% %   'location','best');
% % xlabel('# of eigenvalues'); ylabel('Eigenvalue of M^{-1}L');
% % title(['First 100 eigenvalues of |Y32| perturbed sphere' ...
% %   'with different mesh refinement']);
% % saveas(gcf,'Y32_spectrum.png');
% 
%% cMCF debug
%   load i2_500_t2_abs(Y33(v))_e0.5p0.512.mat
%   numeig = ceil(.5*numv);
  load i3_300_t3_bunny2k_e20p0.mat
  numv = size(v,1); % number of vertices
  numf = size(f,1); % number of faces
  numeig = ceil(.2*numv);
  
%   [M_T,L_T] = lapbel(v_T,f_T);
%   D_T = eigvf(L_T,M_T,numeig);
%   [s_T,v0] = meancurvflow(v_T,f_T,1,'c');
  v0 = v;
  %% prepare to flow it back
  % when is there an edge (mild redundancy)
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

  % compute initial edge lengths squared
  elsq0 = zeros(numv);
  temp = mod(isedge-1,numv)+1;
  temp2 = fix((isedge-1)/numv);
  for j = 1:numv
    for i = temp(temp2 == (j-1))'
      elsq0(i,j) = sum((v(i,:)-v(j,:)).^2);
    end
  end
  elsq0 = nonzeros(elsq0);
  %% re-embedding
  conf_T = sqrt(kron(1./s_T',1./s_T));
  elsq_T = elsq0.*conf_T(isedge); % linear indices
  test = @(v) conformalcost(v,isedge,elsq_T);
  options = optimset('GradObj','on','display','iter-detailed',...
    'maxiter',imax,'largescale','off','tolfun',eps,'tolx',eps);
  [v_Thist,Jc_hist] = fminunc(test,reshape(v',[],1),options);
% [Jc_Thist,v_Thist] = gradescent(@conformalcost,imax,aC,bC,tC,etolC,0,...
%   reshape(v0',[],1),isedge,elsq_T);
  v_c = reshape(v_Thist(:,end),3,[])';

  %% how well does it do?
  [M_c,L_c] = lapbel(v_c,f);
  D_c = eigvf(L_c,M_c,numeig);
  
  viewlim = [-1 1];
  close all;
  figure(); hold all; set(gcf,'outerposition',[0, 0, 1024, 768]);
  subplot(2,3,1); hold all; view(3); grid on; axis equal
  trimesh(f,v_T(:,1),v_T(:,2),v_T(:,3))
  set(gca,'xlim',viewlim,'ylim',viewlim,'zlim',viewlim);
  xlabel('x'); ylabel('y'); zlabel('z');
  title('initial mesh');
  view([-20 80])

  subplot(2,3,2); hold all; view(3); grid on; axis equal
  trisurf(f,v0(:,1),v0(:,2),v0(:,3),s_T,'edgecolor','none')
  set(gca,'xlim',viewlim,'ylim',viewlim,'zlim',viewlim);
  xlabel('x'); ylabel('y'); zlabel('z');
  % colorbar;
  title('S^{-1}');
  text(0,2.5,2,num2str(std(vnorm(v0)),'std(|v|) = %g'));
  view([-20 80])
  
  subplot(2,3,3); hold all; view(3); grid on; axis equal
  trimesh(f,v_c(:,1),v_c(:,2),v_c(:,3))
  set(gca,'xlim',viewlim,'ylim',viewlim,'zlim',viewlim);
  xlabel('x'); ylabel('y'); zlabel('z');
  title('re-embedded mesh');
  view([-20 80])
  
  subplot(2,3,4:6); hold all; grid on;
  v = (D_c - D_T)./D_T;
  plot(v(1:end-1),'ro:','linewidth',2);
  legend('(\lambda_{cMCF embed} - \lambda_{0})/\lambda_{0}',...
    'location','northwest');
  xlabel('# of eigenvalues (#1 is of the highest frequency)');
  title('Deviation from target Laplacian eigenvalues');

  fig = gcf;
  if isstruct(fig)
    ax = fig.CurrentAxes;
  else
    ax = gca;
  end
  pause(.1);
  set(ax,'xlim',[0 numel(D_0)])
  xm = get(ax,'xlim');
  ym = get(ax,'ylim');
  text(floor(max(xm)/3),max(ym(2) + .18*diff(ym)),...
    num2str([numel(Jc_Thist) Jc_Thist(end)],'iter#%d: J_{re-embed} = %g'));
  colormap('jet');

  save('cMCF test.mat','v0','f','v_T','f_T','v_c','D_c',...
    's_T','D_T','Jc_hist');
  hgexport(gcf,'cMCF re-embed test.png',...
    hgexport('factorystyle'), 'Format', 'png'); 

%% line search testing/debug
% close all;
% imax = 1e3;
% % J = @(x) deal(-abs(x),-(x>0)+(x<0));
% J = @(x) deal(abs(x),(x>0)-(x<0));
% % J = @(x) deal(abs(x).^.3,((x>0)-(x<0))*abs(x)^(-2/3)/3);
% % J = @(x) deal(x.^2,2*x);
% % J=@(x)deal(4*x.^3-(x-2).^2+.4*(x+2).^4,12*x.^2-2*(x-2)+1.6*(x+2).^3);
% x0 = -118;
% [Jhist,vhist] = gradescent(J,imax,.5,.6,1,1e-8,0,x0);
% % [Jhist,vhist] = gradescent(J,imax,.7,.8,150,1e-8,0,x0);
% figure(); hold all; grid on;
% xx = [linspace(-abs(x0),0,1e5) linspace(0,abs(x0),1e5)];
% [Jxx,dJxx] = J(xx);
% plot(xx,Jxx,'-');
% plot(vhist,Jhist,'x-','markersize',20)

% conclusion: the smoother the function, the lower c2 needs to be
%% matlab sparse eigs performance test
% ssize = [300 500 1000];
% numeig = .1:.1:.9;
% rat = zeros(numel(ssize),numel(numeig));
% for i = 1:numel(ssize)
%   load([num2str(ssize(i)) '.mat']);
%   for j = 1:numel(numeig)
%     [M,L] = lapbel(v,f);
%     tic
%     [V,lambda,ix] = eigvf(L,M,ceil(numeig(j)*ssize(i)));
%     t1 = toc;
%     M = sparse(M);
%     L = sparse(L);
%     tic
%     [V,lambda,ix] = eigvf(L,M,ceil(numeig(j)*ssize(i)));
%     t2 = toc;
%     rat(i,j) = t2/t1;
%   end
% end
% 
% [X,Y] = meshgrid(ssize,numeig);
% figure(); hold all; grid on; view([60 35]);
% surf(X,Y,rat');
% xlabel('# vtx'); ylabel('amount of eigenvalues requested');
% zlabel('ratio of performance time');