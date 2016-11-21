clear; close all;
%% load target and mcf data and get lapbel
% target = 'spot';
% target = 'spot10k';
% target = 'bunny';
% target = 'spotSCALED';
% target = 'bunnySCALED';
% target = 'spotS10';
% target = 'sphere_large';
% target = 'sphere21k';
% target = 'blob1k';
% target = 'blob4k';
target = 'blob18k';
% target = 'blob96k';

% discriptor = ''; %this is wrong
discriptor = 'inv'; % do not change this
% discriptor = 'log'; %this is also wrong

maxLcomputed = 30; % maxL precomputed in precomputeRSHI.m

Lrange = 6;%[1:3 5 8 10 15 20 25 maxLcomputed];

if exist(['mcf/' target '.mat'],'file')
  load(['mcf/' target '.mat'])
else
  fid = fopen(['../meshes/' target '.obj'],'rt');
  if fid == -1 % try .mat formats
    load(['../meshes/' target '.mat']);
  else
    [v_T,f_T] = readwfobj(fid);
  end

  [s_T,v,sa_T] = meancurvflow(v_T,f_T,1e4,'c');
  save(['mcf/' target '.mat'],'v_T','f_T','s_T','v','sa_T');
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

numv = size(v_T,1);
numf = size(f_T,1);
%% center mesh and normalize surface area

center_v_T = volCenter(v_T,f_T);
v_T = v_T - repmat(center_v_T,numv,1);
v_T = v_T*makeUnitArea(v_T,f_T)*sqrt(4*pi);

center_v = volCenter(v,f_T);
v = v - repmat(center_v,numv,1);
v = v*makeUnitArea(v,f_T)*sqrt(4*pi);

%%

sv = zeros(numv,1);
sa = zeros(numf,1);
area_v_T = zeros(numf,1);
area_v = zeros(numf,1);
for fi = 1:size(f_T,1)
  vi = f_T(fi,:);
  area_v_T(fi) = norm(cross(v_T(vi(2),:) - v_T(vi(1),:), v_T(vi(3),:) - v_T(vi(1),:)))/2;
  area_v(fi) = norm(cross(v(vi(2),:) - v(vi(1),:), v(vi(3),:) - v(vi(1),:)))/2;
  sa(fi) = (area_v(fi)/area_v_T(fi));
end

[M_T,L_T] = lapbel(v_T,f_T);
[M,L] = lapbel(v,f_T);

for fi = 1:size(f_T,1)
  for idx = 0:2
    i = f_T(fi,idx+1);
    j = f_T(fi,mod(idx+1,3)+1);
    k = f_T(fi,mod(idx+2,3)+1);
    sv(i) = sv(i) + sa(fi)*norm(cross(v(j,:) - v(i,:), v(k,:) - v(i,:)))/6;
  end
end
sv = sv./diag(M);

% s_T = sv;

% if strcmp('inv',discriptor)
% elseif strcmp('log',discriptor)
%   s_T = log(1./s_T);
% else
%   s_T = 1./s_T; % re-invert / restore conf. factors
% end

[vham,fham] = xyz2hammer(v);
%%
apjlist = [];
alslist = [];
for idx = 1:numel(Lrange)
maxL = Lrange(idx);
%% get spherical harmonics on vertices
[Y_v,lm] = sphericalHarmonicBase(v,maxL);
numSH = size(lm,1);

%% calculate conf. factors in SH basis
x = s_T;

v_bary = zeros(numf,3); % for 1 point barycentric quadrature
for fi = 1:size(f_T,1)
  vi = f_T(fi,:);
  v_bary(fi,:) = sum(v(vi,:),1)/3;
end

a_pj = zeros((maxL+1)^2,1);
Yvb = sphericalHarmonicBase(v_bary,maxL); % 1 point quadrature
for fi = 1:size(f_T,1)
  a_pj = a_pj + sa(fi)*Yvb(fi,:)'*area_v(fi); % 1 point quadrature
%   a_pj = a_pj + sa(fi)*sum(Y_v(f_T(fi,:),:),1)'*area_v(fi)/3; % PL quadrature
end

a_ls = (Y_v'*Y_v)\(Y_v'*x);
s_ls = Y_v*a_ls;

% delta = Y_v'*M*Y_v;
% a_pj = delta\(Y_v'*M*x);
s_pj = Y_v*a_pj;


if min(s_pj)<=0
  warning('has negative or infinite conformal factors')
  
  az = 0:0.1:2*pi;
  el = 0:0.1:pi;
  [AZ,EL] = meshgrid(az,el);
  R = ones(size(AZ));
  [vx,vy,vz] = sph2cart(AZ(:),EL(:)-pi/2,R(:));
  vv = [vx, vy, vz];
  Y = sphericalHarmonicBase(vv,maxL);
  shsh = Y*a_pj;

  min(shsh)
end

%%
% if strcmp('log',discriptor)
%   modifier = @(x) exp(x);
% elseif strcmp('inv',discriptor)
  modifier = @(x) log(abs(1./x))*(numel(find(x<0,1)) == 1)+ log((1./x))*isempty(find(x<0,1));
% else
%   modifier = @(x) x;
% end

crange = [min(modifier(s_pj)) max(modifier(s_pj))];
figure(); colormap jet
set(gcf,'outerposition',[0, 0, 1920, 1080]);

subplot(2,2,1); hold all; axis equal; view([0 90]);
title('cMCF conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_T),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('eastoutside');
ylabel(ch,'log(1/s)');

subplot(2,2,2); hold all; axis equal; view([0 90]);
title('cMCF mesh grid (hammered)');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3));
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('westoutside');
ylabel(ch,'log(1/s)');

subplot(2,2,3); hold all; axis equal; view([0 90]);
title('SH projected conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_pj),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('eastoutside');
ylabel(ch,'log(1/s)');

subplot(2,2,4); hold all; axis equal; view([0 90]);
title('SH least-squares fitted conformal factors');
trisurf(fham,vham(:,1),vham(:,2),vham(:,3),modifier(s_ls),'facecolor','interp','edgecolor','none');
set(gca,'visible','off'); set(findall(gca, 'type', 'text'), 'visible', 'on');

caxis(crange); ch = colorbar('westoutside');
ylabel(ch,'log(1/s)');
%%
hgexport(gcf,num2str(maxL,['SH/forward/' target 'SphericalHarmonicsL=%02d' discriptor '.png']),...
  hgexport('factorystyle'), 'Format', 'png');
pause(.1);
apjlist = padcat(apjlist,a_pj);
alslist = padcat(alslist,a_ls);

%% test if spherical harmonics look like spherical harmonics 
%  (this debugging procedure failed me because of all the symmetries of SH)
% a = zeros((maxL+1)^2,1);
% a(5) = 2*sqrt(pi);
% sh = Y_v*a;
% 
% figure();
% trisurf(f_T,v(:,1),v(:,2),v(:,3),sh,'facecolor','interp','edgecolor','none'); axis equal;
% colorbar; colormap jet

%% get original spectrum
numeig = 961;
D_T = eigvf(L_T,M_T,numeig);
D_c = eigvf(L,diag(sparse(1./s_T))*M,numeig);

%% get SH basis spectrum
numL = (maxL+1)^2;

D_s = [];
for l = 0:maxLcomputed
    D_s = [D_s, repmat(-l*(l+1),1,2*l+1)]; % spherical harmonic eigenvalues
end

for ll = maxLcomputed:-2:maxL
% for ll = maxLcomputed
numeig = (ll+1)^2;
% numeig = (maxLcomputed+1)^2;
a_pj = [apjlist(1:numL,idx); zeros(numeig-numL,1)];
a_ls = [alslist(1:numL,idx); zeros(numeig-numL,1)];

L_pj = zeros(numeig);
L_ls = zeros(numeig);
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
      apjkcijk = 0;
      alskcijk = 0;
      for k = 1:numeig
        cijk = integrate3realSH(lm([i j k],:));
        apjkcijk = apjkcijk + a_pj(k)*cijk;
        alskcijk = alskcijk + a_ls(k)*cijk;
      end
      L_pj(i,j) = D_s(i)*apjkcijk;
      L_ls(i,j) = D_s(i)*alskcijk;
    end
  else
    L_pj(i,:) = D_s(i)*a_pj'*real(cijk(1:numeig,1:numeig));
    L_ls(i,:) = D_s(i)*a_ls'*real(cijk(1:numeig,1:numeig));
  end
end

D_pj = sort(eig(L_pj));
D_ls = sort(eig(L_ls));

if ll==maxLcomputed
  figure(); hold all; grid on; set(gca,'yscale','log');
end

% plot(numL:-1:1,-D_pj(end-numL+1:end)+D_T(end-numL+1:end))
plot(numeig-1:-1:1,-(-D_pj(end-numeig+1:end-1)+D_T(end-numeig+1:end-1))./D_T(end-numeig+1:end-1))
% set(gca,'yscale','log');
end
% set(gca,'yscale','linear');
legend(num2str((maxLcomputed:-2:maxL)'),'location','best');
title(num2str(maxL,'maxL=%d'));
%% compare spectrum
rei = numL:-1:2;
figure();hold all;
plot(rei,-D_pj(end-numL+1:end-1),'--','linewidth',2)
plot(rei,-D_ls(end-numL+1:end-1),':','linewidth',2)

plot(rei,-D_T(end-numL+1:end-1),'--','linewidth',2)
plot(rei,-D_c(end-numL+1:end-1),'k:','linewidth',2)
% plot(fliplr(rei),-D_s(2:end),':','linewidth',4)
% if strfind(target,'spot')
%   load spotspec.mat
% elseif strfind(target,'bunny')
%   load bunnyspec.mat
% else
%   warning('no discrete spectrum data available');
%   D = [];
% end
% if size(D,1) > numel(rei)
%   rei = [nan(size(size(D,1):-1:numL+1)) rei];
% else
%   rei = rei(end-size(D,1)+1:end-1);
% end
% for i = 1:size(D,2)
%   plot(rei,-D(1:end-1,i));
% end
% set(gca,'xscale','log','yscale','log');
ylim([0 max(-D_pj(end-numL+1:end-1))*1.1])
xlabel('# of Eigenvalues');
ylabel('Laplacian Eigenvalues');
title(num2str(maxL,'Conformal Factor Forward Problem in Spherical Harmonics L=%02d'));
legend('projected SH spectrum','least squares SH spectrum',...
  'Original PL spectrum','McMCF PL spectrum',...
  'location','best');
% saveas(gcf,'SH/forward/spec.png');
saveas(gcf,num2str(maxL,['SH/forward/' target 'SHspecL=%d.png']));
%%
% close all;
figure(); hold all; grid on;
plot(numL:-1:1,abs((-D_pj(end-numL+1:end)+D_T)./D_T));
set(gca,'yscale','log');
xlabel('# of Eigenvalues');
ylabel('Relative Error');
saveas(gcf,num2str(maxL,['SH/forward/' target 'SHspecLogDiffL=%d.png']));

save(num2str(maxL,['SH/forward/' target 'SHspecL=%d.mat']),'a_pj','a_ls','D_pj','D_ls','D_T');
end
% unix(['rm '% save(['SH/forward/' target 'SphericalHarmonics.' discriptor '.mat'],'apjlist','alslist');

%% add caption and make gif when multiple maxL
if numel(Lrange) > 2
close all;
unix(['mogrify -font Liberation-Sans -fill white -undercolor ''#000000F0'' -pointsize 26 ' ...
  '-gravity NorthEast -annotate +10+10 %t SH/forward/' target '*phericalH*' discriptor '.png']);

unix(['convert -delay 100 -loop 0 SH/forward/' target '*phericalH*' ...
  discriptor '.png SH/forward/' target 'SphericalHarmonics.' discriptor '.gif']);

unix(['rm SH/forward/' target '*phericalH*' discriptor '.png']);
end
%% sphere-only test
if strfind(target,'sphere')
figure(); hold all; grid on;
plot(flipud(D_pj)-D_s')
plot(flipud(D_ls)-D_s')
xlabel('# of Eigenvalues');
ylabel('Absolute Error');
title('Sphere spectral error (SH basis - analytic)')
legend('projected','least squares','location','best');
saveas(gcf,num2str(maxL,['SH/forward/SphereSHspecDiffL=%d.png']));
end