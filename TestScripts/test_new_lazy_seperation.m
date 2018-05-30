%first test DB 1
clear all
close all
dbstop if error
dataset = 1;
tic
dt = string(Configuration.TorsoDatasets);
dtt = dt(dataset);
[initial_model, X_ct, minimum_bone_intensity] = SegmentBoneV2(char(dtt), 'end_after_fuzzy_find_bone', 1, 'fixed_initial_threshold',1175);
[~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(char(dtt), false);

dil_map = imdilate(initial_model,[1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1]);

CC = bwconncomp(dil_map);

labelled_initial_model = labelmatrix(CC);
labelled_dilated_model = labelled_initial_model;
labelled_initial_model(initial_model==0) = 0;
%ViewSeperateVolumes(labelled_initial_model,X_ct);

%throw away small CCs
indices_CC = 1:length(CC.PixelIdxList);
for i_CC = 1:length(CC.PixelIdxList)
    if(length(CC.PixelIdxList{i_CC}) < 1000)
        labelled_initial_model(CC.PixelIdxList{i_CC}) = 0;
        labelled_dilated_model(CC.PixelIdxList{i_CC}) = 0;
        indices_CC(indices_CC==i_CC) = [];
    end
end

if(dataset==1)
    bed_inds = [1, 6, 828, 785, 15];
    for i_bedRemove = 1:length(bed_inds)
        labelled_initial_model(CC.PixelIdxList{bed_inds(i_bedRemove)}) = 0;
        labelled_dilated_model(CC.PixelIdxList{bed_inds(i_bedRemove)}) = 0;
        indices_CC(indices_CC==bed_inds(i_bedRemove)) = [];
    end
end

new_model = labelled_initial_model;
new_model(new_model~=0) = 1;
new_CC = bwconncomp(new_model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate edge map
dog_gf1 = 0.0104;
dog_gf2 = 1.9997;
dog_filt = imgaussfilt3(double(ct3D),dog_gf1)-imgaussfilt3(double(ct3D),dog_gf2);
tic
smoothed = imgaussfilt3(dog_filt,4,'FilterSize',11);
toc

dog_minus_smoothed = dog_filt - smoothed;

basic_surfaces = dog_minus_smoothed>-50;

basic_surfaces(ct3D<950) = 0;

bone_plus_surfaces = new_model;
bone_plus_surfaces(basic_surfaces==1) = 1;
bone_plus_surfaces(ct3D<minimum_bone_intensity) = 0;

%fill internal points
%find internal bone
summed_internal_bone = zeros(size(initial_model));
input_mask = imdilate(bone_plus_surfaces,  [1 1 1 ; 1 1 1 ; 1 1 1 ]);

for i_internalPoits = 1:5
    tested_points = input_mask;
    [~,internal_map, n_dirs] = InternalPointsInline(input_mask,bone_plus_surfaces,30,1);
    summed_internal_bone(internal_map==1) = 1;
    
    input_mask = imdilate(internal_map, [1 1 1 ; 1 1 1 ; 1 1 1 ]);
    input_mask(tested_points==1) = 0;
    
end
input_mask = imdilate(summed_internal_bone, [1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1]);
input_mask(summed_internal_bone==1) = 0;
[~,internal_map, num_dir_bone_found_in] = InternalPointsInline(input_mask,bone_plus_surfaces,10,2);
summed_internal_bone(internal_map==1) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%
updated_map = new_model;
updated_map(summed_internal_bone==1) = 1;
map_coating = imdilate(updated_map,[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1]);
map_coating(updated_map==1) = 0;



%map_plus_coating = updated_map;
%map_plus_coating = uint8(map_plus_coating);
%map_plus_coating(map_coating==1) = 1;


%seperate map plus coating
seperated_map = updated_map;
seperated_map(basic_surfaces==1) = 0;

CC_seperated = bwconncomp(seperated_map);
labelled_seperated=labelmatrix(CC_seperated);
% figure
% lengths = [];
% for i_CC = 1:length(CC_seperated.PixelIdxList)
%     lengths(i_CC) = length(CC_seperated.PixelIdxList{i_CC});
% end
% histogram(lengths, 1:length(CC_seperated.PixelIdxList));
updated_labelled = labelled_seperated;
for i_CC = 1:length(CC_seperated.PixelIdxList)
    if(length(CC_seperated.PixelIdxList{i_CC}) < 100)
        updated_labelled(CC_seperated.PixelIdxList{i_CC}) = 0;
        
    end
end

% image = updated_labelled;
% figure
% imshow3D(image,[0 max(max(max(image)))]);


coating_label = max(updated_labelled(:))+1;
updated_labelled(map_coating==1) = coating_label;
fully_labelled = updated_labelled;


%KNN tagging 
%re-fill new labelled map
    labelled_map = zeros(size(updated_labelled));
    labelled_map = uint8(labelled_map);

    labelled_map(updated_labelled~=0) = 1;
    unlabelled_map = uint8(updated_map)-labelled_map;


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
    lookup_n = updated_labelled(labelled_point_indices);

    fully_labelled(unlabelled_point_indices) = lookup_n(IDX);


%remove those tagged as coating
fully_labelled(fully_labelled==coating_label) = 0;

%classify!
labels = unique(fully_labelled); %51 labels!
binary_map = fully_labelled;
binary_map(binary_map>0) = 1;

%GOOD HERE
locality = strel('sphere',2); 
    loc_filter_results = imfilter(binary_map, double(locality.Neighborhood));
    %filter_sum = sum(sum(sum(locality.Neighborhood)));
    %mf = 0.8;
    %rem_mask_ind = loc_filter_results > mf*filter_sum;
    rem_mask = zeros(size(ct3D));
    rem_mask = uint8(rem_mask);
    %rem_mask = uint8(rem_mask);
    %rem_mask(rem_mask_ind) = 1;
    %SE = strel('sphere', size(locality.Neighborhood,1));
    %rem_mask = imdilate(rem_mask, SE);

    %NEW I-1
    rem_mask(binary_map==1) = 1;


    perf_regions = loc_filter_results;
    perf_regions(rem_mask == 1) = 0;
    perf_regions(perf_regions<5) = 0;
    perf_regions(perf_regions~= 0) = 1;
    perf_regions = imdilate(perf_regions, strel('sphere',7));
    perf_regions = uint8(perf_regions);
    


add_bone_edge_map = dog_filt>40;
add_bone_edge_map = uint8(add_bone_edge_map);
add_bone_edge_map(perf_regions==0) = 0;
add_bone_edge_map(ct3D<1100) = 0;

upd_binary_map = binary_map;
upd_binary_map(add_bone_edge_map==1) = 1;

CC_upd = bwconncomp(upd_binary_map);
labelled_m = labelmatrix(CC_upd);
for i_cc = 1:length(CC_upd.PixelIdxList)
    if(length(CC_upd.PixelIdxList{i_cc})<100)
        upd_binary_map(CC_upd.PixelIdxList{i_cc}) = 0;
    end
end



%fill internal points
fix_regions_mask_wo_bone = perf_regions;
fix_regions_mask_wo_bone(upd_binary_map==1) = 0;
[upd_binary_map,internal_map] = InternalPointsInline(fix_regions_mask_wo_bone,upd_binary_map, 10, 1);

for i = 1:5
    fix_regions_mask_wo_bone = imdilate(internal_map, [1 1 1 ; 1 1 1 ; 1 1 1 ]);
    fix_regions_mask_wo_bone(upd_binary_map==1) = 0;
    [upd_binary_map,internal_map] = InternalPointsInline(fix_regions_mask_wo_bone,upd_binary_map, 10, 1);

end

for i = 1:6
    fix_regions_mask_wo_bone = imdilate(upd_binary_map, [1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1]);
    fix_regions_mask_wo_bone(upd_binary_map==1) = 0;
    [upd_binary_map,internal_map] = InternalPointsInline(fix_regions_mask_wo_bone,upd_binary_map, 5, 2);

end

binary_map_filled = imfill(upd_binary_map,'holes');

smoothed = imgaussfilt3(double(binary_map_filled),0.4);

smoothed_thresh = smoothed>0.8;
CC_final = bwconncomp(smoothed_thresh);
for i_cc = 1:length(CC_final.PixelIdxList)
    if(length(CC_final.PixelIdxList{i_cc})<500)
        binary_map_filled(CC_final.PixelIdxList{i_cc}) = 0;
    end
end

smoothed = imgaussfilt3(double(binary_map_filled),0.4);


iso = isosurface(X_ct{1}, X_ct{2}, X_ct{3},smoothed,0.95);
fp = strcat('C:\Users\mazna\Documents\nl\U\P\Data\TestingSets\Torso\8\','bone_map',datestr(now, 'HH-MM-dd-mmm-yyyy'),'.ply');
plywrite(fp,iso.faces,iso.vertices);
