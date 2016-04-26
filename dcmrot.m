function v_result = dcmrot(v_old,dcm)
v_result = v_old';
for i = 1:size(v_old,1)
  v_result(:,i) = dcm*v_old(i,:)';
end
v_result = v_result';