function [hpat] = ViewSeperateVolumes(labelled_map, X_ct, x1_range, x2_range, x3_range)
%VIEWSEPERATEVOLUMES 
%labelled_map fully labelled map (i,e, 1 for first cv, 2 for 2nd and so on
view_map = labelled_map;
%view_map(view_map==0) = nan;
num_cv = max(labelled_map(:));
n_colours = 6;
cmap = jet(n_colours);
for i = 1:num_cv
    full_cmap(i,:) = cmap(mod(i-1,n_colours)+1,:);
end
if(nargin>2)
    [hpat] = PATCH_3Darray(view_map,X_ct{1}(x1_range),X_ct{2}(x2_range),X_ct{3}(x3_range),full_cmap,'col');


else
    [hpat] = PATCH_3Darray(view_map,X_ct{1},X_ct{2},X_ct{3},full_cmap,'col');


end

end

