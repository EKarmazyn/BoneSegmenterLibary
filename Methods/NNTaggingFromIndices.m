function [taggedIdx_array] = NNTaggingFromIndices(seedIdx_array, to_tag_indices, X_ct, img_size)
%NNTAGGING 

all_seed_indexes = [];
lookup_n =[];
for i_seed = 1:length(seedIdx_array)
    all_seed_indexes = [all_seed_indexes; seedIdx_array{i_seed}];
    lookup_n = [lookup_n, ones(1,length(seedIdx_array{i_seed}))*i_seed];
end

%fully_labelled = labelled_seed_map;

%KNN tagging 
%re-fill new labelled map
    %labelled_map = zeros(size(to_tag_map));
    %labelled_map = uint8(labelled_map);

    %labelled_map(all_seed_indexes) = 1;
    %unlabelled_map = uint8(to_tag_map);
    %unlabelled_map(labelled_map==1) = 0;


    %find nearest point in labelled map for each point in
    %unlabelled_map
    unlabelled_point_indices = to_tag_indices;
    
    %remove tagged indices
    unlabelled_point_indices = setdiff(unlabelled_point_indices,all_seed_indexes);
    
    [ulp_1,ulp_2,ulp_3] = ind2sub(img_size,unlabelled_point_indices);

    %labelled_point_indices = find(labelled_map);
    labelled_point_indices = all_seed_indexes;
    [lp_1,lp_2,lp_3] = ind2sub(img_size,labelled_point_indices);
    %convert to real positions using x_ct
    ulp_1 = X_ct{1}(ulp_1);
    ulp_2 = X_ct{2}(ulp_2);
    ulp_3 = X_ct{3}(ulp_3);

    lp_1 = X_ct{1}(lp_1);
    lp_2 = X_ct{2}(lp_2);
    lp_3 = X_ct{3}(lp_3);

    [IDX, D] = knnsearch([lp_1',lp_2',lp_3'],[ulp_1',ulp_2',ulp_3']);

    
    
    
    %IDX maps unlabelled_point_indices to labelled_point_indices, which are
    %in the seedIdx array
    
    % create lookup for matching
    %lookup_n = labelled_seed_map(labelled_point_indices);
    mapping = lookup_n(IDX);
    
    taggedIdx_array = seedIdx_array;
    for i_seed = 1:length(seedIdx_array)
        taggedIdx_array{i_seed} = [seedIdx_array{i_seed}; unlabelled_point_indices(mapping==i_seed)];
    end
    
    
    %fully_labelled(unlabelled_point_indices) = lookup_n(IDX);
    %fully_labelled = reshape(fully_labelled,size(labelled_seed_map));

end

