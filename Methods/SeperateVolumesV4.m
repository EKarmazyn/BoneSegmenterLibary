function [new_connected_volumes_list, new_labelled_map] = SeperateVolumesV4(connected_volumes,ct3D, X_ct, dog_filt, knn_filling, dog_filt_threshold)
%SEPERATEVOLUMES4 utilises a mesh -> hole filling -> binary approach
addpath(genpath('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\PotentiallyUsefulLibraries\iso2mesh'));
tic
new_connected_volumes_list = [];

new_labelled_map = zeros(size(dog_filt));
new_labelled_map = uint16(new_labelled_map);

voided_map = zeros(size(ct3D));
voided_map = uint8(voided_map);
voided_map(dog_filt>dog_filt_threshold) = 1;

%hole filling!
%flat surface 
    
for i_cc = 1:length(connected_volumes)
    
    bin_map = zeros(size(ct3D));
    bin_map = uint8(bin_map);
    bin_map(connected_volumes(i_cc).IndVoxelList) = 1;
    
    
    
    %masked_dog = zeros(size(dog_filt));
    %masked_dog = uint8(masked_dog);
    %masked_dog(connected_volumes(i_cc).IndVoxelList) = voided_map(connected_volumes(i_cc).IndVoxelList);
    
   
    %
    bin_map2 = bin_map;
    bin_map2(voided_map==1) = 0;
    
    [CC] = bwconncomp(bin_map2);
    %[labelled_binary_mapCC,N] = bwlabeln(bin_map2);
    labelled_binary_mapCC = labelmatrix(CC);
    
    if(isempty(CC.PixelIdxList))
        new_connected_volumes_list = [new_connected_volumes_list; connected_volumes(i_cc)];
        new_labelled_map(connected_volumes(i_cc).IndVoxelList) = length(new_connected_volumes_list);
        continue;
    end
    
    
    if(knn_filling)
        %re-fill new labelled map
        labelled_map = zeros(size(ct3D));
        labelled_map = uint8(labelled_map);
        labelled_map(labelled_binary_mapCC==1) = 1;
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
        lookup_n = labelled_binary_mapCC(labelled_point_indices);
        
        labelled_binary_mapCC(unlabelled_point_indices) = lookup_n(IDX);
        
    end
    for n = 1:length(CC.PixelIdxList)
        
        new_connected_volumes_list = [new_connected_volumes_list; ConnectedVolume(find(labelled_binary_mapCC==n))];
        new_labelled_map(labelled_binary_mapCC==n) = length(new_connected_volumes_list);
    end
end






toc



end

