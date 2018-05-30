function [connected_volumes_array, fully_labelled] = PreClassificationSeperate(ct3D,binary_map_filled, X_ct, CT_dimension_spacing, xy_thresh,z_thresh)
%PRECLASSIFICATIONSEPERATE Summary of this function goes here
%seperate step---------------------------------------


voided = HardSeperateVolumes(ct3D,binary_map_filled, CT_dimension_spacing, xy_thresh,z_thresh);
%[connected_volumes2, labelled_binary_map2,all_points_ind_list2] = ConnectedVolumes(voided);

%looks good

CC = bwconncomp(voided);
labelled_seed_map = labelmatrix(CC);

min_hard_bone = 1500;

%KNN tag the HARD BONE first (since some is lost entirely (i.e. bits of rib
%/ ilium
hard_bone_mask = zeros(size(ct3D));
hard_bone_mask(ct3D>min_hard_bone) = 1;
hard_bone_mask(binary_map_filled==0) = 0;

[hard_labelled] = NNTagging(labelled_seed_map, hard_bone_mask, X_ct);


%coating
map_coating = imdilate(binary_map_filled,[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1]);
map_coating(binary_map_filled==1) = 0;


coating_label = max(hard_labelled(:))+1;
hard_labelled(map_coating==1) = coating_label;

%KNN tag rest
[fully_labelled] = NNTagging(hard_labelled, binary_map_filled, X_ct);

%remove coating
fully_labelled(fully_labelled==coating_label) = 0;
labels = unique(fully_labelled); %18328 labels!

updated_binary_map = fully_labelled;
updated_binary_map(updated_binary_map>0) = 1;


labelled_indices = find(updated_binary_map);
indices_labels = fully_labelled(labelled_indices);

label_index_list = cell(length(labels),1);

for i_label = 1:length(labels)
     %tic
        label_index_list{i_label,1} = labelled_indices(indices_labels==labels(i_label));
     %toc
end

valid_labelled_index_list = cell(length(label_index_list)-1,1);

for i_label = 1:length(labels)-1
     %tic
        valid_labelled_index_list{i_label,1} = label_index_list{i_label+1};
     %toc
end

connected_volumes_array = ConnectedVolume.FromLabelledIndices(valid_labelled_index_list);



end

