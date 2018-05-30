function [connected_volumes, labelled_binary_map,all_points_ind_list] = ConnectedVolumes(binary_map)
%CONNECTEDVOLUMES Generates a list of connected volumes from a binary map



[labelled_binary_map,N] = bwlabeln(binary_map);

CC = bwconncomp(binary_map);
reorder_array = nan(N,1); %%maps CC COs to bwlabeln labels (for re-ordering the connected volume array before function end)

for co = 1:N
    
    %cur_co = L(L==co);
    %check lengths match
    
    cur_co = CC.PixelIdxList{co}; 
    obj_length = length(cur_co);
    connected_volumes(co) = ConnectedVolume(cur_co);
    
    %find the corresponding label and put in in the reorder array
    reorder_array(co) = labelled_binary_map(cur_co(1));
    
    
end
[~,invSort] = sort(reorder_array);
connected_volumes = connected_volumes(invSort);



%CC = bwconncomp(binary_map);
%connected_volumes = nan(CC.NumObjects,1);
% for co = 1:CC.NumObjects
%     cur_co = CC.PixelIdxList{co};
%     obj_length = length(cur_co);
%     connected_volumes(co) = ConnectedVolume(cur_co);
% end
temp_cell_array = {};
for j = 1:length(CC.PixelIdxList)
    temp_cell_array{j,1} = CC.PixelIdxList{j};
end
all_points_ind_list = cell2mat(temp_cell_array);

end

