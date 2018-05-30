function [new_labelled_map] = SeperateVolumesV6(connected_volumes,ct3D,cur_bone_model_map, X_ct, dog_gf1, dog_gf2, dog_filt_threshold, seperation_use_model_surface,  seperation_min_voxels, dilate_join_internals)
%SEPERATEVOLUMES3

%tic
%new_connected_volumes_list = [];

%r = 1;
new_labelled_map = zeros(size(ct3D));
new_labelled_map = uint16(new_labelled_map);


dog_filt = imgaussfilt3(double(ct3D),dog_gf1)-imgaussfilt3(double(ct3D),dog_gf2);


surfaces_map = dog_filt>dog_filt_threshold;
surfaces_map = uint8(surfaces_map);
surfaces_map(cur_bone_model_map==0) = 0;
CC = bwconncomp(surfaces_map);
for i_cc = 1:length(CC.PixelIdxList)
    if(length(CC.PixelIdxList{i_cc})<seperation_min_voxels)
        surfaces_map(CC.PixelIdxList{i_cc}) = 0;
    end
end

for i_slice = 1:size(ct3D,3)
    %tic
    filled_slice = imfill(surfaces_map(:,:,i_slice),'holes');
    dil = imdilate(surfaces_map(:,:,i_slice),[1 1 1; 1 1 1]);
    dil(filled_slice==1) = 0;
    surfaces_map(:,:,i_slice) = surfaces_map(:,:,i_slice) + dil;
    %round 2
    filled_slice = imfill(surfaces_map(:,:,i_slice),'holes');
    dil = imdilate(surfaces_map(:,:,i_slice),[1 1 1; 1 1 1]);
    dil(filled_slice==1) = 0;
    surfaces_map(:,:,i_slice) = surfaces_map(:,:,i_slice) + dil;
    %toc
end

%add edge of bone model to surfaces map?
if(seperation_use_model_surface)
    surfaces_map = surfaces_map + (uint8(cur_bone_model_map) - uint8(imerode(cur_bone_model_map,[1 1 1; 1 1 1; 1 1 1])));
    surfaces_map(surfaces_map==2) = 1;
end
    
new_map_counter = 1;

for i_cc = 1:length(connected_volumes)
    
    bin_map = zeros(size(ct3D));
    bin_map = uint8(bin_map);
    bin_map(connected_volumes(i_cc).IndVoxelList) = 1;
    bin_map2 = bin_map;
    bin_map2(surfaces_map==1) = 0;
    
    %[base_bound_map,ind1,ind2,ind3, full_size] = BoundMap(bin_map,2);
    %[bound_map,ind1,ind2,ind3, full_size] = BoundMap(bin_map2,2);
    bound_map = bin_map2;
    if(isempty(bound_map))
        %new_connected_volumes_list = [new_connected_volumes_list; connected_volumes(i_cc)];
        new_labelled_map(connected_volumes(i_cc).IndVoxelList) = new_map_counter;
        new_map_counter = new_map_counter+1;
        continue;
    end
    
    %re-dilate based off of surface thickness
    if(dilate_join_internals)
        bound_map = imdilate(bound_map, [1 1 1; 1 1 1; 1 1 1]);
    end
    
    [CC] = bwconncomp(bound_map);
    %[labelled_binary_mapCC,N] = bwlabeln(bin_map2);
    labelled_binary_mapCC = uint16(labelmatrix(CC));
    
    if(length(CC.PixelIdxList)==0)
        new_labelled_map(connected_volumes(i_cc).IndVoxelList) = new_map_counter;
        new_map_counter = new_map_counter+1;
        continue;
    end
    
    %full_size_labelledCC = UnBoundMap(labelled_binary_mapCC, ind1,ind2,ind3, full_size);
    full_size_labelledCC = labelled_binary_mapCC;
    
    
    %re-fill new labelled map
    labelled_map = zeros(size(full_size_labelledCC));
    labelled_map = uint8(labelled_map);

    labelled_map(full_size_labelledCC~=0) = 1;
    unlabelled_map = bin_map-labelled_map;


    %find nearest point in labelled map for each point in
    %unlabelled_map
    unlabelled_point_indices = find(unlabelled_map);
    [ulp_1,ulp_2,ulp_3] = ind2sub(size(unlabelled_map),unlabelled_point_indices);

    labelled_point_indices = find(labelled_map);
    [lp_1,lp_2,lp_3] = ind2sub(size(unlabelled_map),labelled_point_indices);
    %convert to real positions using x_ct
    ulp_1 = X_ct{1}(ulp_1);
    ulp_2 = X_ct{2}(ulp_2);
    ulp_3 = X_ct{3}(ulp_3);

    lp_1 = X_ct{1}(lp_1);
    lp_2 = X_ct{2}(lp_2);
    lp_3 = X_ct{3}(lp_3);

    [IDX, D] = knnsearch([lp_1',lp_2',lp_3'],[ulp_1',ulp_2',ulp_3']);


    % create lookup for matching
    lookup_n = full_size_labelledCC(labelled_point_indices);

    full_size_labelledCC(unlabelled_point_indices) = lookup_n(IDX);
        
    
    
    
        
    %new_connected_volumes_list = [new_connected_volumes_list; ConnectedVolume(find(full_size_labelledCC==n))];
    new_labelled_map(full_size_labelledCC~=0) = full_size_labelledCC(full_size_labelledCC~=0)+new_map_counter-1;
    new_map_counter = new_map_counter+length(CC.PixelIdxList);
     
end






%toc



end

