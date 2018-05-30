function [] = SegmentBoneV1(varargin)
%SEGMENTBONEV1 

data_filepath = varargin{1,1};
if(length(varargin) > 1)
    options = varargin(2:length(varargin));
end


%% Set Default options
seperate_radius = 1;
seperate_erode = true;



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
[Slices, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(data_filepath, false);


%% Valid Area Mask
[valid_mask] = CreateValidAreaMask(ct3D);

%% Initial Straight Threshold
[bone_model_map, minimum_bone_intensity] = FuzzyFindBone(ct3D, valid_mask,true, false);

%% Perforated Fix

locality = strel('sphere',2); 
loc_filter_results = imfilter(bone_model_map, double(locality.Neighborhood));
%filter_sum = sum(sum(sum(locality.Neighborhood)));
%mf = 0.8;
%rem_mask_ind = loc_filter_results > mf*filter_sum;
rem_mask = zeros(size(ct3D));
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

gray_fix_map = ct3D;
gray_fix_map(perf_regions~=1) = 0;
gray_fix_map_wo_bone = gray_fix_map;
gray_fix_map_wo_bone(bone_model_map==1) = 0;

gf1 = imgaussfilt3(ct3D,0.5);
gf2 = imgaussfilt3(ct3D,2);
dog_filt_full = gf1-gf2;
dog_filt = dog_filt_full;
dog_filt(perf_regions==0) = 0;
dog_filt(bone_model_map==1) = 0;

add_bone_edge_map = zeros(size(gray_fix_map));
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

%% Fill inernal points 1

fix_regions_mask_wo_bone = perf_regions;
fix_regions_mask_wo_bone(cur_model==1) = 0;
[updated_bone_map,internal_map] = InternalPoints(fix_regions_mask_wo_bone,cur_model, true, 20);

%NEW I-3
%% Fill Internal Points 2 

%% CV Analysis
[connected_volumes, labelled_binary_map, all_points_ind_list] = ConnectedVolumes(updated_bone_map);

%% Sub-seperate volumes
[connected_volumes, new_labelled_map] = SeperateVolumesV3(connected_volumes, ct3D, X_ct, dog_filt_full, 1, 0, true,20);


%% Classify

%find probable and definate bone seeds
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


SaveBone3D(new_labelled_map, [bone_vols; not_bone_vols], X_ct);

%% PSS to mesh



end

