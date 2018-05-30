function [new_labelled_map, X_ct, minimum_bone_intensity] = SegmentBoneV3(varargin)
%SEGMENTBONEV3, like V2 but internal-filling after seperation

data_filepath = varargin{1,1};
if(length(varargin) > 1)
    options = varargin(2:length(varargin));
end


%% Set Default options
seperation_dilation_passes = 2;
dog_gf1 = 0.0104;
dog_gf2 = 1.9997;
seperation_min_voxels = 6;
dog_filt_threshold = 0;
seperation_use_model_surface = true;
end_after_seperation = false;
end_after_fuzzy_find_bone = false;
no_debug_messages = false;
perf_regions = true; %limit perf_fix to perf_identified regions


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
[valid_mask] = CreateValidAreaMask(ct3D);
if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

%% Initial Straight Threshold
if(~no_debug_messages)
    fprintf("FuzzyFindBone");
    tic
end
[bone_model_map, minimum_bone_intensity] = FuzzyFindBone(ct3D, valid_mask);
close;
if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end
fprintf("fuzzy_find_bone_intensity: " + num2str(minimum_bone_intensity));
if(end_after_fuzzy_find_bone)
    new_labelled_map = [];
    return;
end
%% Perforated Fix

if(~no_debug_messages)
    fprintf("Perforated Fix");
    tic
end

gray_fix_map = ct3D;

if(perf_regions)
    locality = strel('sphere',2); 
    loc_filter_results = imfilter(bone_model_map, double(locality.Neighborhood));
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
    rem_mask(bone_model_map==1) = 1;


    perf_regions = loc_filter_results;
    perf_regions(rem_mask == 1) = 0;
    perf_regions(perf_regions<5) = 0;
    perf_regions(perf_regions~= 0) = 1;
    perf_regions = imdilate(perf_regions, strel('sphere',7));
    perf_regions = uint8(perf_regions);

    %gray_fix_map = ct3D;
    gray_fix_map(perf_regions~=1) = 0;
    %gray_fix_map_wo_bone = gray_fix_map;
    %gray_fix_map_wo_bone(bone_model_map==1) = 0;
else
    perf_regions = uint8(ones(size(ct3D)));
end
gf1 = imgaussfilt3(ct3D,0.5);
dog_filt = gf1-imgaussfilt3(ct3D,2);
dog_filt(perf_regions==0) = 0;
dog_filt(bone_model_map==1) = 0;

add_bone_edge_map = zeros(size(gray_fix_map));
add_bone_edge_map = uint8(add_bone_edge_map);
dog_thresh = 20;
add_bone_edge_map(dog_filt>dog_thresh) = 1;

%NEW I-2
add_bone_edge_map(ct3D<970) = 0;
add_bone_edge_map(gf1<1120) = 0;

cur_model = add_bone_edge_map+bone_model_map;
cur_model(cur_model==2) = 1;

CC = bwconncomp(cur_model);
for i = 1:length(CC.PixelIdxList)
    if(length(CC.PixelIdxList{i}) < 1000)
        cur_model(CC.PixelIdxList{i}) = 0;
        
    end
end

if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

%% CV Analysis
if(~no_debug_messages)
    fprintf("CV Analysis");
    tic
end
[connected_volumes, labelled_binary_map, all_points_ind_list] = ConnectedVolumes(cur_model);
if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

%% Erode Pre-Seperate chaf

%% Sub-seperate volumes
if(~no_debug_messages)
    fprintf("Sub-seperate volumes");
    tic
end

[new_labelled_map] = SeperateVolumesV5(connected_volumes,ct3D, cur_model, X_ct, dog_gf1,dog_gf2, dog_filt_threshold, seperation_use_model_surface, seperation_min_voxels);
if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end



%% CV Analysis 2
if(~no_debug_messages)
    fprintf("CV Analysis 2");
    tic
end
connected_volumes = [];
for i_cv = 1:max(new_labelled_map(:))
    connected_volumes = [connected_volumes; ConnectedVolume(new_labelled_map==i_cv)];
end

if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end


%% Fill inernal points 1 (on seperate connected volumes version)


if(~no_debug_messages)
    fprintf("Fill inernal points 1");
    tic
end
cur_cv_mask = uint8(zeros(size(ct3D)));
for i_cv = 1:length(connected_volumes)
    cur_cv(connected_volumes(i_cv).IndVoxelList) = 1;
    cur_cv_mask = CustomDilate(cur_cv,10);
    [~,internal_map] = InternalPoints(cur_cv_mask,cur_cv, true, 20);
    new_labelled_map(internal_map==1) = i_cv;
end

if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end


%% CV analysis 3
if(~no_debug_messages)
    fprintf("CV Analysis 3");
    tic
end
connected_volumes = [];
for i_cv = 1:max(new_labelled_map(:))
    connected_volumes = [connected_volumes; ConnectedVolume(new_labelled_map==i_cv)];
end
if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

if(end_after_seperation)
    return
end

%NEW I-3
%% Fill Internal Points 2 


%% Classify

%find probable and definate bone seeds
if(~no_debug_messages)
    fprintf("Classify");
    tic
end
[prob_bone_seed_map, prob_minimum_bone_seed_intensity, def_bone_seed_map, def_bone_seed_min_intensity] = FindBoneSeeds(ct3D, valid_mask);

%mark connected volumes with bone seeds
ConnectedVolume.MarkWithBoneSeeds(connected_volumes,prob_bone_seed_map, def_bone_seed_map);



%extract features of connected volumes
for cv = 1:length(connected_volumes)
    connected_volumes(cv).GenerateIndividualFeatures(ct3D);
end


ConnectedVolume.GenerateRegionPropsFeatures(connected_volumes, new_labelled_map);


%load model
load('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\Models\PC3\new_classi_model2GOOD.mat');

%temp_ind = [1:11, 13:length(connected_volumes)];
[full_table, bone_table, non_bone_table] = ConnectedVolume.ConstructDataTable(connected_volumes);
%input_table = [bone_table; non_bone_table];

%test model
yfit = new2_classi_model.predictFcn(full_table);
tested_connected_volumes = connected_volumes;
%Hightlight3D(fuzzy_bone_map,connected_volumesALT(3),[],[]);
bone_ind = (yfit==1);
not_bone_ind = (yfit==-1);
%Highlight3D(fuzzy_bone_map,tested_connected_volumes(not_bone_ind),tested_connected_volumes(bone_ind),[]);


bone_vols = tested_connected_volumes(bone_ind);
for bv = 1:length(bone_vols)
    bone_vols(bv).Bone = 1;
end

not_bone_vols = tested_connected_volumes(not_bone_ind);
for bv = 1:length(not_bone_vols)
    not_bone_vols(bv).Bone = -1;
end
if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

if(~no_debug_messages)
    fprintf("Saving");
    tic
end
SaveBone3D(new_labelled_map, [bone_vols; not_bone_vols], X_ct);
if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end
%% PSS to mesh



end

