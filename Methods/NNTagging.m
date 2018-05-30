function [fully_labelled] = NNTagging(labelled_seed_map, to_tag_map, X_ct)
%NNTAGGING 

fully_labelled = labelled_seed_map;

%KNN tagging 
%re-fill new labelled map
    labelled_map = zeros(size(labelled_seed_map));
    labelled_map = uint8(labelled_map);

    labelled_map(labelled_seed_map~=0) = 1;
    unlabelled_map = uint8(to_tag_map);
    unlabelled_map(labelled_map==1) = 0;


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
    lookup_n = labelled_seed_map(labelled_point_indices);

    fully_labelled(unlabelled_point_indices) = lookup_n(IDX);
    %fully_labelled = reshape(fully_labelled,size(labelled_seed_map));

end

