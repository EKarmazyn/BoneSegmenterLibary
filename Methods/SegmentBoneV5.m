function [connected_volumes_array, X_ct, minimum_bone_intensity] = SegmentBoneV5(varargin)
%SEGMENTBONEV4 new_lazy_seperation


data_filepath = varargin{1,1};
if(length(varargin) > 1)
    options = varargin(2:length(varargin));
end


%% Set Default options
xy_thresh = 15;
z_thresh = 0;
pre_volume_classifier_save = false;
seperation_dilation_passes = 2;
dog_gf1 = 0.0104;
dog_gf2 = 1.9997;
seperation_min_voxels = 6;
dog_filt_threshold = 0;
seperation_use_model_surface = false;
end_after_seperation = false;
end_after_fuzzy_find_bone = false;
end_after_filling = false;
no_debug_messages = false;
perf_regions = true; %limit perf_fix to perf_identified regions
dilate_join_internals = true;
added_bone_min_voxels = 100;
fixed_initial_threshold = 0;
dataset = -1;
save_results = false;


%% Process Options for overrides

if(nargin > 1)
    num_options = length(options)/2;

    for i = 1:num_options
        o_i = 2*i-1;
        feval(@()assignin('caller',options{o_i}, options{o_i+1}));
        %assignin('caller',options{o_i},options{o_i+1});
        
    end
end





%% Load DICOM (intensity data)
if(~no_debug_messages)
    fprintf("LoadingDICOM");
    tic
end
[~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(data_filepath, false);
if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end


%% Valid Area Mask
if(~no_debug_messages)
    fprintf("CreateValidAreaMask");
    tic
end
%[valid_mask] = CreateValidAreaMask(ct3D);
valid_mask = ones(size(ct3D));
valid_mask = uint8(valid_mask);
if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

%% Initial Straight Threshold
if(~no_debug_messages)
    fprintf("Generate Initial Model");
    tic
end
if(fixed_initial_threshold==0)
    [bone_model_map, minimum_bone_intensity] = FuzzyFindBone(ct3D, valid_mask);
    close;
else
    bone_model_map = uint8(ct3D>fixed_initial_threshold);
    minimum_bone_intensity = fixed_initial_threshold;
end

%fprintf("fuzzy_find_bone_intensity: " + num2str(minimum_bone_intensity));
if(end_after_fuzzy_find_bone)
    new_labelled_map = bone_model_map;
    return;
end

%bone_model_map
initial_model = bone_model_map;


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

new_model = labelled_initial_model;
new_model(new_model~=0) = 1;
new_CC = bwconncomp(new_model);

if(save_results)
    SaveResults(new_model, "initial_model_" + num2str(dataset));
end

if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

%% SS1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~no_debug_messages)
    fprintf("SS1");
    tic
end

%generate edge map
dog_gf1 = 0.0104;
dog_gf2 = 1.9997;
dog_filt = imgaussfilt3(double(ct3D),dog_gf1)-imgaussfilt3(double(ct3D),dog_gf2);
%tic
smoothed = imgaussfilt3(dog_filt,4,'FilterSize',11);
%toc

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



if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SS2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~no_debug_messages)
    fprintf("SS2");
    tic
end

%
updated_map = new_model;
updated_map(summed_internal_bone==1) = 1;

if(save_results)
    SaveResults(updated_map, "post_ss1_" + num2str(dataset));
end

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


if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(~no_debug_messages)
    fprintf("SS3");
    tic
end


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

if(save_results)
    SaveResults(fully_labelled, "post_ss2_" + num2str(dataset));
end


if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SS3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(~no_debug_messages)
    fprintf("SS4");
    tic
end



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

if(save_results)
    SaveResults(upd_binary_map, "post_ss3_" + num2str(dataset));
end


%% SS4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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

%binary_map_filled = imfill(upd_binary_map,'holes');


if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(pre_volume_classifier_save)
    %save .mat binary map
    fprintf("Saving... \n");
    save("C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\PreSeperation\" + num2str(dataset),'binary_map_filled')
    
    return;
end

if(save_results)
    SaveResults(upd_binary_map, "post_ss4_" + num2str(dataset));
end



%% Volume classifier

%HardSeperate----------------------------------------

%hole filling? YES!
eroded_bin_map = imerode(upd_binary_map, CrossFilter(3));

CC_eroded = bwconncomp(eroded_bin_map);

all_ind = find(upd_binary_map);
to_upd_ind = setdiff(find(eroded_bin_map), all_ind);

%tag unlabelled
[taggedIdx_array] = NNTaggingFromIndices(CC_eroded.PixelIdxList,to_upd_ind, X_ct, size(ct3D));

zero_map = uint16(zeros(size(ct3D)));
binary_map_filled = uint8(zero_map);

for i_cv = 1:length(taggedIdx_array)
    cur_map = zero_map;
    cur_map(taggedIdx_array{i_cv}) = 1;
    binary_map_filled(find(imfill(cur_map,'holes'))) = 1;
    
end
 %imfill(upd_binary_map,'holes');

%no scale:
SCALE;
SCALE.xy = 0.792968750000000;
SCALE.z = 1.500000000000000;
SCALE.mean = (1.500000000000000 + 2*0.79296875)/3;
SCALE.vol = 0.943199157714844;


[labelled_map] = PreClassifySeperate2(ct3D, binary_map_filled, CT_dimension_spacing, X_ct,SCALE, false);


pre_process_lab_map = labelled_map;
[connected_volumes_array] = PostPreClassifySeparateProcess(pre_process_lab_map, ct3D,dataset,X_ct,SCALE);

if(save_results)
    SaveResults(connected_volumes_array, "pre_classify_seperate" + num2str(dataset));
end


%Classify----------------------------------------




%% Clean

%% Convert

%% Save


end

