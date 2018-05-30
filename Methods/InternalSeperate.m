function [indVoxelLists] = InternalSeperate(connected_volume, full_size, max_points)
%INTERNALSEPERATE Splits a connected volume into a set of smaller regions
%connected volume: the CV to sub-split
%full size - size of ct3D
%max_points - number of points to aim for with each internally seperated
%region
%indVoxelLists - returned as indices of the full map not the bound (as borders differ) region!
indVoxelLists = {};


cv_map_full = zeros(full_size);
cv_map_full(connected_volume.IndVoxelList) = 1;


[cv_map_bound,ind1,ind2,ind3, ~] = BoundMap(cv_map_full,1);
inv_map = find(cv_map_bound == 0);

unused_voxel_list = connected_volume.IndVoxelList;
unused_region = cv_map_bound;




start_point =unused_voxel_list(1);
[sp_sub1,sp_sub2,sp_sub3] = ind2sub(full_size,start_point);
sp_bound_sub = [sp_sub1-(ind1-1), sp_sub2-(ind2-1), sp_sub3-(ind3-1)];
cur_region = zeros(size(cv_map_bound));
sp_bound = sub2ind(size(cv_map_bound), sp_bound_sub(1),sp_bound_sub(2),sp_bound_sub(3));
cur_region(sp_bound)=1;


i = 1;
el = 15;
while(true)
    
    %start point:
    if(i~=1)
        start_point =unused_voxel_list(1);
        cur_region = zeros(size(cv_map_bound));
        cur_region(start_point)=1;
    end
    SE = strel('cube',el);
    region_length = 1;
    while(region_length < max_points)
        prev_length = region_length;
        %can boundary on current region
        [cr_map_bound,crind1,crind2,crind3, cr_size_full] = BoundMap(cur_region,el);
        
        cr_map_bound = imdilate(cr_map_bound, SE);
        
        [cur_region] = UnBoundMap(cr_map_bound,crind1,crind2,crind3,cr_size_full);
        
        %mask
        cur_region(inv_map) = 0;
        region_length = sum(sum(sum(cur_region)));

        if(region_length == prev_length)
            break;
        end
    end
    indVoxelLists{i} = find(cur_region==1);
    
    %update unused_voxel_list:
    unused_region(find(cur_region==1)) = 0;
    unused_voxel_list = find(unused_region);
    toc
    tic
    if(length(unused_voxel_list) == 0)
        break;
    end
    inv_map = find(unused_region==0);
    
    i = i+1;
    
end

%for regions with length less than half max join them to other regions?


%convert voxelLists to full map indices
for i_vl = 1:length(indVoxelLists)
    cur_in_bound = indVoxelLists{i_vl};
    bs_map = zeros(size(cv_map_bound));
    bs_map(cur_in_bound) = 1;
    on_full = UnBoundMap(bs_map,ind1,ind2,ind3,full_size);
    indVoxelLists{i_vl} = find(on_full);
end


end

