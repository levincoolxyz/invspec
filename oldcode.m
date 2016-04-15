% regular octahedron (1)
if input_case == 1
  v = [0 0 1; 0 1 0; 1 0 0 ;-1 0 0; 0 -1 0; 0 0 -1];
  f=fliplr(convhulln(v));
  
% regular icosahedron (2)
% cf. Anton Semechko (a.semechko@gmail.com)
elseif input_case == 2
  t=(1+sqrt(5))/2; % golden ratio
  v=[0 1 t];
  s=[1 1 1; 1 1 -1; 1 -1 -1; 1 -1 1];
  v=repmat(v,[4 1]).*s;
  v=[v;circshift(v,[0 -1]);circshift(v,[0 -2])];
  v_L2=sqrt(sum(v.^2,2));
  v=bsxfun(@rdivide,v,v_L2);
  clear s; clear t;
  f=fliplr(convhulln(v));
end

% skewed octahedron
if target_case == 1
  v_T = [0 0 .5; 0 1 0; 2 0 0 ;-2 0 0; 0 -3 0; 0 0 -.5];
  f_T = fliplr(convhulln(v));
  [M_T,L_T] = lapbel(v_T,f_T);
  D_T = eigvf(L_T,M_T,numeig);
end