function scl = makeUnitArea(v,f)
area = 0;
for fi = 1:size(f,1)
  vi = f(fi,:);
  area = area + norm(cross(v(vi(2),:) - v(vi(1),:),...
    v(vi(3),:) - v(vi(1),:)))/2;
end
scl = 1./sqrt(area);
end