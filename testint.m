% syms x;
% assume(x,'real');
% %%
% IC = [];
% for c = -5:.1:5
%   IC = [IC;eval(int(sin(x)*sin(c*x),x,0,2*pi))];
% end
%%
% plot(0:.1:5,IC)

x = -5:.1:5;

plot(x,(abs(x)~=1).*sin(2*pi*x)./(x.^2-1)+pi*(abs(x)==1))