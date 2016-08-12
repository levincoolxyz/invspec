clear;
%% load target and mcf data
% target = 'spot';
% target = 'spot10k';
target = 'bunny';
% target = 'spotS10';
% target = 'sphere_large';

% discriptor = '';
discriptor = 'inv';
% discriptor = 'log';

if exist(['mcf/' target '.mat'],'file')
  load(['mcf/' target '.mat'])
else
  fid = fopen(['../meshes/' target '.obj'],'rt');
  if fid == -1 % try .mat formats
    load(['../meshes/' target '.mat']);
  else
    [v_T,f_T] = readwfobj(fid);
  end
  [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');
  save(['mcf/' target '.mat'],'v_T','f_T','s_T','v');
end

% pert = @(x) .5*(2*v(:,3).^2-v(:,1).^2-v(:,2).^2)./(vnorm(v)).^2;
% 
% [M,L] = lapbel(v_T,f_T);
% invM = diag(1./diag(M));
% Hn = .5*[invM*L*v(:,1) invM*L*v(:,2) invM*L*v(:,3)];
% H = vnorm(Hn);
% vn = Hn./repmat(H,1,3);
% v_T = v_T - repmat(pert(v),1,3).*vn;
% [s_T,v] = meancurvflow(v_T,f_T,1e5,'c');

if strcmp('inv',discriptor)
elseif strcmp('log',discriptor)
  s_T = log(1./s_T);
else
  s_T = 1./s_T; % re-invert / restore conf. factors
end

[vham,fham] = xyz2hammer(v);
%%
apjlist = [];
alslist = [];
for maxL = 10%[1:3 5 8 10 15 20 25 40]
%% get spherical harmonics on vertices
[Y_v,LM] = sphericalHarmonicBase(v,maxL);
numSH = size(LM,1);

%% calculate conf. factors in SH basis
x = s_T;

a_ls = (Y_v'*Y_v)\(Y_v'*x);
s_ls = Y_v*a_ls;

M = lapbel(v,f_T);
delta = Y_v'*M*Y_v;
a_pj = delta\(Y_v'*M*x);
s_pj = Y_v*a_pj;

%%
if strcmp('log',discriptor)
  modifier = @(x) exp(x);
elseif strcmp('inv',discriptor)
  modifier = @(x) log(abs(1./x))*(numel(find(x<0,1)) == 1)+ log((1./x))*isempty(find(x<0,1));
else
  modifier = @(x) x;
end

crange = [min(modifier(s_pj)) max(modifier(s_pj))];
figure(); colormap jet
set(gcf,'outerposition',[0, 0, 1920, 1080]);

subplot(2,2,1); hold all; axis equal; view([0 90]);
title('cMCF conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_T),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('eastoutside');
ylabel(ch,'s');

subplot(2,2,2); hold all; axis equal; view([0 90]);
title('cMCF mesh grid (hammered)');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3));
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('westoutside');
ylabel(ch,'s');

subplot(2,2,3); hold all; axis equal; view([0 90]);
title('SH projected conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_pj),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('eastoutside');
ylabel(ch,'s');

subplot(2,2,4); hold all; axis equal; view([0 90]);
title('SH least-squares fitted conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_ls),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('westoutside');
ylabel(ch,'s');
%%
% hgexport(gcf,num2str(maxL,[target 'SphericalHarmonicsL=%02d' discriptor '.png']),...
%   hgexport('factorystyle'), 'Format', 'png');
pause(.1);
apjlist = padcat(apjlist,a_pj);
alslist = padcat(alslist,a_ls);
end
% save(['SH/' target 'SphericalHarmonics.' discriptor '.mat'],'apjlist','alslist');

%% add caption and make gif
% unix(['mogrify -font Liberation-Sans -fill white -undercolor ''#000000F0'' -pointsize 26 ' ...
%   '-gravity NorthEast -annotate +10+10 %t ' target '*phericalH*' discriptor '.png']);
% 
% unix(['convert -delay 100 -loop 0 ' ...
%   target '*phericalH*' discriptor '.png SH/' target 'SphericalHarmonics.' discriptor '.gif']);
% 
% unix(['rm ' target '*phericalH*' discriptor '.png']);

%% test if spherical harmonics look like spherical harmonics
% a = zeros((maxL+1)^2,1);
% a(5) = 2*sqrt(pi);
% sh = Y_v*a;
% min(sh)
% % min(s_pj)
% 
% % az = 0:0.1:2*pi;
% % el = 0:0.1:pi;
% % [AZ,EL] = meshgrid(az,el);
% % R = ones(size(AZ));
% % [vx,vy,vz] = sph2cart(AZ(:),EL(:)-pi/2,R(:));
% % vv = [vx, vy, vz];
% % Y = sphericalHarmonicBase(vv,10);
% % % shsh = Y*a;
% % shsh = Y*a_pj;
% % % shsh = zeros(size(vx));
% % % for j = 1:numel(vx)
% % %   for i = 1:121
% % %     shsh(j) = shsh(j) + Y(j,i)*a_pj(i);
% % %   end
% % % end
% % 
% % min(shsh)
% 
% figure();
% trisurf(f_T,v(:,1),v(:,2),v(:,3),sh,'facecolor','interp','edgecolor','none'); axis equal;
% colorbar; colormap jet

%% get original spectrum
numeig = max(maxL+1)^2;
[M_T,L_T] = lapbel(v_T,f_T);
D_T = eigvf(L_T,M_T,numeig);

% [M,L] = lapbel(v,f_T);
% D_w = eigvf(L,diag(sparse(s_T)).*M,numeig);

%% get SH basis spectrum
% a = zeros(size(a_pj));
% a(1) = 2*sqrt(pi);
% aa = a;
% aa = [2*sqrt(pi), -1.3182409, 0.66036477, -0.34151179, 0.20749373, -0.083409802 ...
% 0.0082154643, -0.025159468, -0.0062353568, -0.0029661895 ...
% 0.0011572657, 0.0011554944, 0.00066818418, 0.00021389462 ... 
% -0.00014639729, 0.000072971848]';
aa = a_pj;
% aa = a_ls;

D_s = [];
for l = 0:maxL
    D_s = [D_s, repmat(-l*(l+1),1,2*l+1)]; % spherical harmonic eigenvalues
end

L_sh = zeros(numeig);
for i = 1:numeig
  if exist(num2str(i,'../RSHI/RSHI%04d.mat'),'file')
    load(num2str(i,'../RSHI/RSHI%04d.mat'));
    computeNow = 0;
  else
    warning('ain''t nobody got time for this');
    computeNow = 1;
  end
  if computeNow
    for j = 1:numeig
      akcijk = 0;
      for k = 1:numeig
        cijk = integrate3realSH(LM([i j k],:));
        akcijk = akcijk + aa(k)*cijk;
      end
      L_sh(i,j) = D_s(i)*akcijk;
    end
  else
    L_sh(i,:) = D_s(i)*aa'*real(cijk(1:numeig,1:numeig));
  end
end

D_sh = eig(L_sh);
% if mean(abs(imag(D_sh))) > 1e-10
%   error('no no no no no no ...');
% end
% D_sh = sort(real(D_sh));
D_sh = sort(D_sh);

%% compare spectrum
rei = numeig:-1:1;
figure();hold all;
plot(rei,-D_sh,'--')
%%
plot(rei,-D_T)
if strfind(target,'spot')
  load spotspec.mat
elseif strfind(target,'bunny')
  load bunnyspec.mat
else
  warning('no discrete spectrum data available');
  D = [];
end
if size(D,1) > numel(rei)
  rei = [nan(size(size(D,1):-1:numeig+1)) rei];
else
  rei = rei(end-size(D,1)+1:end);
end
for i = 1:size(D,2)
  plot(rei,-D(:,i));
end