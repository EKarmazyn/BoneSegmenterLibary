function [new_connected_volumes_list, new_labelled_map] = ErodeSeperate(bone_model_map, pre_erode_labelled_map, X_ct, r)
%ERODESEPERATE Re-seperates by simple erosion then knn filling

new_connected_volumes_list = [];

%r = 1;
SE = strel('sphere',r);


eroded_map = imerode(bone_model_map,SE);
[CC] = bwconncomp(eroded_map);
labelled_binary_mapCC = labelmatrix(CC);


%use new labelled_binary_mapCC to update pre_erode_labelled_map copy
%(mid_laballed_map)
mid_laballed_map_old_labels = pre_erode_labelled_map;
mid_laballed_map_old_labels(labelled_binary_mapCC==0) = 0;
mid_laballed_map = zeros(size(mid_laballed_map_old_labels));
zero_map = zeros(size(mid_laballed_map_old_labels));
zero_map = uint16(zero_map);
counter_old_labels = 1;
counter_new_labels = 1;
for i = 1:length(CC.PixelIdxList)
    temp_map = zero_map;
    temp_map(mid_laballed_map_old_labels==counter_old_labels) = 1;
    
    [CC2] = bwconncomp(temp_map);
    temp_labelled_map = labelmatrix(CC2);
    if(max(temp_labelled_map(:))>1)
        %write each into mid_labelled_map with new_labels
        for i2 = 1:length(CC2.PixelIdxList)
           
            mid_laballed_map(temp_labelled_map==i2) = counter_new_labels;
            counter_new_labels = counter_new_labels+1;
        end
        
        %knn filling blanks
        currently_labelled = zero_map;
        

        currently_labelled(temp_labelled_map~=0) = 1;
        unlabelled_map = temp_map-currently_labelled;


        %find nearest point in labelled map for each point in
        %unlabelled_map
        unlabelled_point_indices = find(unlabelled_map);
        [ulp_1,ulp_2,ulp_3] = ind2sub(size(unlabelled_map),unlabelled_point_indices);

        labelled_point_indices = find(currently_labelled);
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
        lookup_n = mid_laballed_map(labelled_point_indices);

        mid_laballed_map(unlabelled_point_indices) = lookup_n(IDX);
        
        %%% KNN FILLING END
        
    else
        mid_laballed_map(mid_laballed_map_old_labels==counter_old_labels) = counter_new_labels;
        counter_new_labels = counter_new_labels+1;
        %no change
    end

    
    
    counter_old_labels = counter_old_labels+1;
end

%knn filling

labelled_map = zeros(size(bone_model_map));
labelled_map = uint8(labelled_map);

labelled_map(labelled_binary_mapCC~=0) = 1; %labelled after erosion
unlabelled_map = bone_model_map-labelled_map;


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
lookup_n = mid_laballed_map(labelled_point_indices);

new_labelled_map = mid_laballed_map;
%new_labelled_map(unlabelled_point_indices) = 0;

new_labelled_map(unlabelled_point_indices) = lookup_n(IDX);

 for n = 1:max(mid_laballed_map(:))
    new_connected_volumes_list = [new_connected_volumes_list; ConnectedVolume(find(new_labelled_map==n))];
    %new_labelled_map(labelled_binary_mapCC==n) = length(new_connected_volumes_list);
end



end

