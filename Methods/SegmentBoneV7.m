function [connected_volumes_array, X_ct] = SegmentBoneV7(varargin)
%SEGMENTBONEV4 uses faster executables for the very slow bits


data_filepath = varargin{1,1};
if(length(varargin) > 1)
    options = varargin(2:length(varargin));
end

%% Set Default options
xy_thresh = 15;
z_thresh = 0;
pre_volume_classifier_save = false;
seperation_dilation_passes = 2;
dog_gf1_mult = 0.0104;
dog_gf2_mult = 1.9997;
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
save_append = "";
ss_save = false;
ss_path = 'C:\Users\mazna\Documents\nl\U\P\Report\FinalReport\FigureData\ss\';
SaveFolder = '';
save_all = false;

sl =300;
kp.x = 120;
kp.y = 120;
kp.z = 160;


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
    bigT = tic;
end
[~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(data_filepath, false);

XY_MM_SCALE = CT_dimension_spacing(1);
Z_MM_SCALE = CT_dimension_spacing(3);
MEAN_MM_SCALE = (2*XY_MM_SCALE + Z_MM_SCALE)/3;

VOXEL_VOLUME = XY_MM_SCALE * XY_MM_SCALE * Z_MM_SCALE;

SCALE.xy = XY_MM_SCALE;
SCALE.z = Z_MM_SCALE;
SCALE.mean = MEAN_MM_SCALE;
SCALE.vol = VOXEL_VOLUME;

%re-scale vars
dog_gf1 = dog_gf1_mult * XY_MM_SCALE * 1.2611;
dog_gf2 = dog_gf2_mult * XY_MM_SCALE * 1.2611;


zero8 = uint8(zeros(size(ct3D)));

if(~no_debug_messages)
    t = toc(bigT);
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end


%% Valid Area Mask
if(~no_debug_messages)
    fprintf("CreateValidAreaMask");
    bigT = tic;
end
%[valid_mask] = CreateValidAreaMask(ct3D);
valid_mask = ones(size(ct3D));
valid_mask = uint8(valid_mask);
if(~no_debug_messages)
    
    t = toc(bigT);
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

%% Initial Straight Threshold
if(~no_debug_messages)
    fprintf("Generate Initial Model");
   bigT = tic;
    
end
if(fixed_initial_threshold==0)
    [bone_model_map, initial_threshold, probable_min_bone_intensity] = FuzzyFindBone(ct3D, valid_mask);
    SCALE.initial_t = initial_threshold;
    close;
else
    bone_model_map = uint8(ct3D>fixed_initial_threshold);
    %minimum_bone_intensity = fixed_initial_threshold;
    initial_threshold = fixed_initial_threshold;
    SCALE.initial_t = fixed_initial_threshold;
end

%fprintf("fuzzy_find_bone_intensity: " + num2str(minimum_bone_intensity));
if(end_after_fuzzy_find_bone)
    new_labelled_map = bone_model_map;
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SAVE
SaveKeyPointsImgs(kp,sl,bone_model_map,SaveFolder+"1init_model")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bone_model_map
initial_model = bone_model_map;

d1_r = round(MEAN_MM_SCALE*1.9443); %2 for dataset 1
dil_map = imdilate(initial_model,CrossFilter(d1_r));

CC = bwconncomp(dil_map);

labelled_initial_model = labelmatrix(CC);
labelled_dilated_model = labelled_initial_model;
labelled_initial_model(initial_model==0) = 0;
%ViewSeperateVolumes(labelled_initial_model,X_ct);

%throw away small CCs
indices_CC = 1:length(CC.PixelIdxList);

min_vol = 2378.9*0.3965;%VOXEL_VOLUME
min_voxels = round(min_vol/VOXEL_VOLUME); %1000 for 1


for i_CC = 1:length(CC.PixelIdxList)
    if(length(CC.PixelIdxList{i_CC}) < (min_voxels/2))
        labelled_initial_model(CC.PixelIdxList{i_CC}) = 0;
        labelled_dilated_model(CC.PixelIdxList{i_CC}) = 0;
        indices_CC(indices_CC==i_CC) = [];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SAVE
SaveKeyPointsImgs(kp,sl,labelled_initial_model,SaveFolder+"2rem_small_vols")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


new_model = labelled_initial_model;
new_model(new_model~=0) = 1;
new_CC = bwconncomp(new_model);

if(save_results)
    SaveResults(new_model, "initial_model_" + save_append + "_" + num2str(dataset));
end

if(~no_debug_messages)
   
    t = toc(bigT);
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

%% SS1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~no_debug_messages)
    fprintf("SS1");
    bigT = tic;
    
end

%generate edge map
dog_gf1_v2 = dog_gf1;
dog_gf2_v2 = dog_gf2;

%note this isnt doing z / xy edges separately because the lower z resultion means they arent clear so shouldnt be found here
dog_filt = imgaussfilt3(double(ct3D),dog_gf1_v2)-imgaussfilt3(double(ct3D),dog_gf2_v2);
%tic


smooth1_filt_size = round(XY_MM_SCALE * 13.8719); %11 for 1
smooth1_filt_size = round(smooth1_filt_size/2);
if(mod(smooth1_filt_size,2)==0)
    smooth1_filt_size = smooth1_filt_size+1;
end
smoothed = imgaussfilt3(dog_filt,4,'FilterSize',smooth1_filt_size);
%toc

dog_minus_smoothed = dog_filt - smoothed;


basic_surfaces = dog_minus_smoothed>-50; %>-50;

%bone surfaces at min will be around 1100
%bone surfaces in the dog map should at min be around 30

dog_minus_smoothed_c = uint16(dog_minus_smoothed) + ct3D;
%basic_surfaces = zero8;
basic_surfaces(dog_minus_smoothed_c<(1100)) = 0;

bone_plus_surfaces = new_model;
bone_plus_surfaces(basic_surfaces==1) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SAVE
SaveKeyPointsImgs(kp,sl,bone_plus_surfaces,SaveFolder+"3ExtraSurfaces")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bone_plus_surfaces(ct3D<probable_min_bone_intensity) = 0;%THIS STEP NEEDS CHECKING

%smooth and threshold?

%level of smoothing can be determined from the amount left over, we are
%looking to add a small amount to the initial model

%init_num = sum(new_model(:));
%extra_num = sum(bone_plus_surfaces(:)) - init_num;
%~1.4 = extra_num/init_num
%we ideally want 0.2 or so extra
%cur_extra(1) = extra_num/init_num;

%smoothed_model = imgaussfilt3(double(bone_plus_surfaces),4,'FilterSize',5);
%bone_plus_surfaces(smoothed_model<0.18) = 0;

smoothed_model = imgaussfilt3(bone_plus_surfaces,4,'FilterSize',7);
bone_plus_surfaces(smoothed_model<0.5) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SAVE
SaveKeyPointsImgs(kp,sl,bone_plus_surfaces,SaveFolder+"4FilteredOut")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

%fill internal points
%find internal bone
summed_internal_bone = zeros(size(initial_model));
dil_r_2 = round(0.9722 * MEAN_MM_SCALE); %1 for 1

input_mask = imdilate(uint8(bone_plus_surfaces),  CrossFilter(dil_r_2));

search_r_1 = round(MEAN_MM_SCALE * 29.166); %30 for 1

comb = bone_plus_surfaces;
summed_internal_bone = bone_plus_surfaces;



for i_internalPoits = 1:4
    tested_points = uint8(input_mask);
    tested_points(summed_internal_bone==1) = 0;
   
    [internal_map] = RemoteInternalFill(tested_points,uint8(summed_internal_bone),search_r_1,5);
    summed_internal_bone(internal_map==1) = 1;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SAVE
    SaveKeyPointsImgs(kp,sl,summed_internal_bone,SaveFolder+"5_" + num2str(i_internalPoits) +"infill"+num2str(i_internalPoits));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
        comb(internal_map==1) = i_internalPoits+1;
    
    input_mask = imdilate(summed_internal_bone, CrossFilter(dil_r_2));
    input_mask(summed_internal_bone==1) = 0;
end

% 
% input_mask = imdilate(summed_internal_bone,CrossFilter(10))-summed_internal_bone;
% 
%   
% for i_internalPoits = 5:9
%     tested_points = uint8(input_mask);
%     tested_points(summed_internal_bone==1) = 0;
%    
%     [internal_map] = RemoteInternalFill(tested_points,uint8(summed_internal_bone),20,13);
%     summed_internal_bone(internal_map==1) = 1;
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%% SAVE
%     SaveKeyPointsImgs(kp,sl,summed_internal_bone,SaveFolder+"5_" + num2str(i_internalPoits) +"infill"+num2str(i_internalPoits));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     
%         comb(internal_map==1) = i_internalPoits+1;
%     
%     input_mask = imdilate(summed_internal_bone, CrossFilter(10));
%     input_mask(summed_internal_bone==1) = 0;
% end
% 
% %these are seeds!
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%% SAVE
%     SaveKeyPointsImgs(kp,sl,comb,SaveFolder+"5COMB_"+"infill"+num2str(i_internalPoits));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     




dil_r_3 = round(MEAN_MM_SCALE * 1.9444); %2 for 1

search_r_2 = round(MEAN_MM_SCALE * 9.722); %10 for 1


input_mask = imdilate(summed_internal_bone, CrossFilter(dil_r_3));
input_mask(summed_internal_bone==1) = 0;


%fprintf("CS_InternalFill");
 %   t =tic;
    [internal_map] = RemoteInternalFill(input_mask,uint8(bone_plus_surfaces),search_r_2,7);
 %   cs_time = toc(t);
%     fprintf(".        Time-Taken: " + num2str(cs_time)+"\n");
   

%[~,internal_map, num_dir_bone_found_in] = InternalPointsInline(input_mask,bone_plus_surfaces,search_r_2,2);
summed_internal_bone(internal_map==1) = 1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SAVE
    SaveKeyPointsImgs(kp,sl,summed_internal_bone,SaveFolder+"6infill_plus1");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(ss_save)
    comb(internal_map==1) = max(comb(:))+1;
    SaveResults(comb, "ss1_2"+ save_append + "_" + num2str(dataset));
end

if(~no_debug_messages)
   
    t = toc(bigT);
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SS2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~no_debug_messages)
    fprintf("SS2");
    bigT = tic;
    
end

%
updated_map = new_model;
updated_map(summed_internal_bone==1) = 1;

if(save_results)
    SaveResults(updated_map, "post_ss1_" + save_append + "_" + num2str(dataset));
end

coating_r_1 = round(2.9166* MEAN_MM_SCALE); %3 for 1

map_coating = imdilate(updated_map,CrossFilter(coating_r_1));
map_coating(updated_map==1) = 0;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SAVE
    SaveKeyPointsImgs(kp,sl,summed_internal_bone,SaveFolder+"7MapCoating");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%map_plus_coating = updated_map;
%map_plus_coating = uint8(map_plus_coating);
%map_plus_coating(map_coating==1) = 1;


%seperate map plus coating
seperated_map = updated_map;
seperated_map(basic_surfaces==1) = 0;



CC_seperated = bwconncomp(seperated_map);
labelled_seperated=labelmatrix(CC_seperated);

if(ss_save)
    
    SaveResults(labelled_seperated, "ss2_1"+ save_append + "_" + num2str(dataset));
end

% figure
% lengths = [];
% for i_CC = 1:length(CC_seperated.PixelIdxList)
%     lengths(i_CC) = length(CC_seperated.PixelIdxList{i_CC});
% end
% histogram(lengths, 1:length(CC_seperated.PixelIdxList));
updated_labelled = labelled_seperated;


min_cc_vol_2 = 237.8906*0.3965;
min_cc_length_2 = round(min_cc_vol_2/VOXEL_VOLUME); %100 for 1

for i_CC = 1:length(CC_seperated.PixelIdxList)
    if(length(CC_seperated.PixelIdxList{i_CC}) < min_cc_length_2)
        updated_labelled(CC_seperated.PixelIdxList{i_CC}) = 0;
        
    end
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SAVE
    SaveKeyPointsImgs(kp,sl,updated_labelled,SaveFolder+"8throwSmall");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(ss_save)
    
    SaveResults(updated_labelled, "ss2_2"+ save_append + "_" + num2str(dataset));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




coating_label = max(updated_labelled(:))+1;
updated_labelled(map_coating==1) = coating_label;
fully_labelled = updated_labelled;


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SAVE
    SaveKeyPointsImgs(kp,sl,fully_labelled,SaveFolder+"9coating");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SAVE
    SaveKeyPointsImgs(kp,sl,fully_labelled,SaveFolder+"10taggedWcoat");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remove those tagged as coating
fully_labelled(fully_labelled==coating_label) = 0;


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SAVE
    SaveKeyPointsImgs(kp,sl,fully_labelled,SaveFolder+"10taggedWOcoat");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%classify!
labels = unique(fully_labelled); %51 labels!
binary_map = fully_labelled;
binary_map(binary_map>0) = 1;

if(save_results)
    SaveResults(fully_labelled, "post_ss2_" + save_append + "_" + num2str(dataset));
end
if(~no_debug_messages)
   
    t = toc(bigT);
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SS3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(~no_debug_messages)
    fprintf("SS3");
   bigT = tic;
    
end




%GOOD HERE
loc_r = round(1.9444*MEAN_MM_SCALE); %2 for 1

locality = strel('sphere',loc_r); 
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
    if(ss_save)
        SaveResults(perf_regions, "ss3_perf_regions"+ save_append + "_" + num2str(dataset));
    end

    dil_3_r = round(6.8054*MEAN_MM_SCALE); %7! for 1
    perf_regions = imdilate(perf_regions, strel('sphere',dil_3_r));
    perf_regions = uint8(perf_regions);
    


add_bone_edge_map = dog_filt>40;
add_bone_edge_map = uint8(add_bone_edge_map);
add_bone_edge_map(perf_regions==0) = 0;
add_bone_edge_map(ct3D<1100) = 0;




upd_binary_map = binary_map;
upd_binary_map(add_bone_edge_map==1) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SAVE
    SaveKeyPointsImgs(kp,sl,upd_binary_map,SaveFolder+"11updBoneEdge");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC_upd = bwconncomp(upd_binary_map);
labelled_m = labelmatrix(CC_upd);
for i_cc = 1:length(CC_upd.PixelIdxList)
    if(length(CC_upd.PixelIdxList{i_cc})<min_cc_length_2)
        upd_binary_map(CC_upd.PixelIdxList{i_cc}) = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SAVE
    SaveKeyPointsImgs(kp,sl,upd_binary_map,SaveFolder+"12throwsmall");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ss_save)
    comb = binary_map;
    comb(add_bone_edge_map==1) = 2;
    comb(upd_binary_map==0) =0;
    SaveResults(comb, "ss3_added_bone"+ save_append + "_" + num2str(dataset));
end
if(save_results)
    SaveResults(upd_binary_map, "post_ss3_" + save_append + "_" + num2str(dataset));
end
if(~no_debug_messages)
   
    t = toc(bigT);
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

%% SS4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~no_debug_messages)
    fprintf("SS4");
    bigT = tic;
   
end


%fill internal points
fix_regions_mask_wo_bone = perf_regions;
fix_regions_mask_wo_bone(upd_binary_map==1) = 0;



fprintf("CS_InternalFill");
    t =tic;
    [internal_map] = RemoteInternalFill(fix_regions_mask_wo_bone,uint8(upd_binary_map),search_r_2,5);
    cs_time = toc(t);
     fprintf(".        Time-Taken: " + num2str(cs_time)+"\n");
   upd_binary_map(internal_map==1) = 1;
%[upd_binary_map,internal_map] = InternalPointsInline(fix_regions_mask_wo_bone,upd_binary_map, search_r_2, 1);


    comb = upd_binary_map;
    comb(internal_map==1) = 2;
    
    


for i = 1:5
    fix_regions_mask_wo_bone = imdilate(internal_map, CrossFilter(dil_r_2));
    fix_regions_mask_wo_bone(upd_binary_map==1) = 0;
    
    
    fprintf("CS_InternalFill");
    t =tic;
    [internal_map] = RemoteInternalFill(fix_regions_mask_wo_bone,uint8(upd_binary_map),search_r_2,5);
    cs_time = toc(t);
     fprintf(".        Time-Taken: " + num2str(cs_time)+"\n");
   upd_binary_map(internal_map==1) = 1;
    
    comb(internal_map==1) = 2+i;
    %[upd_binary_map,internal_map] = InternalPointsInline(fix_regions_mask_wo_bone,upd_binary_map, search_r_2, 1);

end


for i = 1:6
    fix_regions_mask_wo_bone = imdilate(upd_binary_map, CrossFilter(dil_r_3));
    fix_regions_mask_wo_bone(upd_binary_map==1) = 0;
    
    

    [internal_map] = RemoteInternalFill(fix_regions_mask_wo_bone,uint8(upd_binary_map),round(search_r_2/2),8);
  
   upd_binary_map(internal_map==1) = 1;
   
   comb(internal_map==1) = 2+5+i;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SAVE
    SaveKeyPointsImgs(kp,sl,comb,SaveFolder+"15Afterinfilling");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~no_debug_messages)
 
    t = toc(bigT);
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Volume classifier
if(~no_debug_messages)
    fprintf("FinalSeparation");
    tic
end

%HardSeperate----------------------------------------

binary_map_filled = upd_binary_map; %no filling version as was too slow,


[labelled_map, taggedIdx_array] = PreClassifySeperate2(ct3D, binary_map_filled, CT_dimension_spacing, X_ct,SCALE, false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SAVE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pre_process_lab_map = labelled_map;
[connected_volumes_array] = PostPreClassifySeparateProcess(pre_process_lab_map,taggedIdx_array, ct3D,dataset,X_ct,SCALE);



if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end

%Classify---and save

if(~no_debug_messages)
    fprintf("VolumeClassifyAndSave");
    tic
end

[final_cva, binary_map_filled] = VolumeClassifier(connected_volumes_array,X_ct, size(ct3D));

if(~no_debug_messages)
    t = toc;
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SAVE
SaveKeyPointsImgs(kp,sl,binary_map_filled,SaveFolder+"20after_classifier.mat")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save("C:\Users\mazna\Documents\nl\U\P\Data\ML_results\PostClass\FINAL_CVA_"+num2str(dataset)+datestr(now, 'HH-MM-dd-mmm-yyyy') +".mat" ,'final_cva');





end

