function [original_map] = UnBoundMap(bound_map,ind1,ind2,ind3,full_size)
%UNBOUNDMAP undoes bound_map
original_map = zeros(full_size);
l1=  size(bound_map,1)-1;
l2 = size(bound_map,2)-1;
l3 = size(bound_map,3)-1;
original_map(ind1:ind1+l1, ind2:ind2+l2, ind3:ind3+l3) = bound_map;
end



         