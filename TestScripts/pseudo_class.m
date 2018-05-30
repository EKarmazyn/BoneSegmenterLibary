

dataset=1;
dt = string(Configuration.TorsoDatasets);
dtt = dt(dataset);
    
[~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(char(dtt), false);

XY_MM_SCALE = CT_dimension_spacing(1);
Z_MM_SCALE = CT_dimension_spacing(3);
MEAN_MM_SCALE = (2*XY_MM_SCALE + Z_MM_SCALE)/3;

VOXEL_VOLUME = XY_MM_SCALE * XY_MM_SCALE * Z_MM_SCALE;

SCALE.xy = XY_MM_SCALE;
SCALE.x = Z_MM_SCALE;
SCALE.mean = MEAN_MM_SCALE;
SCALE.vol = VOXEL_VOLUME;
    
[connected_volumes_array] = PostPreClassifySeparateProcess(final_labelled_map,ct3D,dataset,X_ct,SCALE );

fully_labelled_map = uint16(zeros(size(ct3D)));
for i = 1:length(connected_volumes_array)
    
     
    
    fully_labelled_map(connected_volumes_array(i).IndVoxelList) = i;
    
        
end


%select group of cvs
% 5 = v_large
vl_labelled_map = uint16(zeros(size(ct3D)));
vl_array = [];
l_labelled_map = uint16(zeros(size(ct3D)));
l_array = [];
m_labelled_map = uint16(zeros(size(ct3D)));
m_array = [];
s_labelled_map = uint16(zeros(size(ct3D)));
s_array = [];
vs_labelled_map = uint16(zeros(size(ct3D)));
vs_array = [];
ss_lab_map = uint16(zeros(size(ct3D)));
cur_size = 5; 

cur_list_counter = 1;
cva_lengths = [];
super_short_list = [];
for i = 1:length(connected_volumes_array)
    cva_lengths = [cva_lengths; connected_volumes_array(i).NumVoxels];
    if(connected_volumes_array(i).VolumeSize == 5)
      
            vl_array = [vl_array; connected_volumes_array(i)];
            vl_labelled_map(connected_volumes_array(i).IndVoxelList) = i;
            %cur_list_counter = cur_list_counter + 1;
      
        
    elseif(connected_volumes_array(i).VolumeSize == 4)
        l_array = [l_array; connected_volumes_array(i)];
            l_labelled_map(connected_volumes_array(i).IndVoxelList) = i;
            %cur_list_counter = cur_list_counter + 1;
    elseif(connected_volumes_array(i).VolumeSize == 3)
        m_array = [m_array; connected_volumes_array(i)];
            m_labelled_map(connected_volumes_array(i).IndVoxelList) = i;
            %cur_list_counter = cur_list_counter + 1;
    elseif(connected_volumes_array(i).VolumeSize == 2)
        s_array = [s_array; connected_volumes_array(i)];
            s_labelled_map(connected_volumes_array(i).IndVoxelList) = i;
            %cur_list_counter = cur_list_counter + 1;
    elseif(connected_volumes_array(i).VolumeSize == 1)
        vs_array = [vs_array; connected_volumes_array(i)];
            vs_labelled_map(connected_volumes_array(i).IndVoxelList) = i;
            %cur_list_counter = cur_list_counter + 1;
    end
    
    if(connected_volumes_array(i).NumVoxels<=5)
        super_short_list = [super_short_list, i];
        connected_volumes_array(i).VolumeSize = -1;
        ss_lab_map(connected_volumes_array(i).IndVoxelList) = i;
    end
end

%purge bad!
valid = 1:length(connected_volumes_array);
valid = setdiff(valid,super_short_list);

full_wo_ss_label = fully_labelled_map;
full_wo_ss_label(ss_lab_map~=0) = 0;


%%FULL ATM
connected_volumes = connected_volumes_array;



%generate features
%mark connected volumes with bone seeds
%[prob_bone_seed_map, prob_minimum_bone_seed_intensity, def_bone_seed_map, def_bone_seed_min_intensity] = FindBoneSeeds(ct3D, ones(size(ct3D)));

%ConnectedVolume.MarkWithBoneSeeds(connected_volumes,prob_bone_seed_map, def_bone_seed_map);


%extract features of connected volumes
for cv = 1:length(connected_volumes)
    connected_volumes(cv).GenerateIndividualFeatures(ct3D);
    
end


ConnectedVolume.GenerateRegionPropsFeatures(connected_volumes, fully_labelled_map);

ConnectedVolume.GenerateSpatialFeatures(connected_volumes, fully_labelled_map,X_ct);


%mark only the worst ones as not bone:
target_vol = super_short_list;
target_vol_inds = find(ss_lab_map~=0);
use_indices = true;

view_individual_volume_script
not_bones = [1288, 91, 139, 171, 1390];
bones = [];
%not_bones = [4 25 292 297 31 267 268 177];

%bones = [331 332 324 326 7 8 11];

%reset all
for i_cv = 1:length(connected_volumes)
    connected_volumes(i_cv).Bone = 0;
end

for i_nb = 1:length(not_bones)
    connected_volumes(not_bones(i_nb)).Bone = -1;
end


for i_b = 1:length(bones)
    connected_volumes(bones(i_b)).Bone = 1;
end

%local solidity features? i.e. how much as % is lost on small erode?
binary_map = full_wo_ss_label;
binary_map(binary_map~=0 ) = 1;

ConnectedVolume.GenerateErodeFraction(connected_volumes_array, binary_map);



[full_table, bone_table, non_bone_table, labelled_table, unlabelled_table] =  ConnectedVolume.ConstructDataTable(connected_volumes);
valid_table = full_table(valid,:);
% figure
% imshow3D(labelled_map);
% 


%colour by erode fraction!

%red-green [0 0 1] -> [1 0 0]
for i_va = 1:length(connected_volumes_array)
    frac = connected_volumes_array(i_va).Features.ErodeFraction;
    frac_a(i_va) = frac;
    size_a(i_va) = connected_volumes_array(i_va).NumVoxels;
    col_arr{i_va} = [ frac (1-frac) 0];
    if(connected_volumes_array(i_va).NumVoxels>800)
        col_arr{i_va} = [ 0 0 1];
        hpat{i_va}.FaceAlpha = 0.1;
        hpat{i_va}.EdgeAlpha = 0.1;
    end
    
    if(connected_volumes_array(i_va).NumVoxels<5)
        col_arr{i_va} = [ 1 0 0];
    end
    
end


erode_valid_ids = find(frac_a<0.7);
final_model = zero_map;
final_valid = intersect(valid,erode_valid_ids);
for i_fv = 1:length(connected_volumes_array)
    if(sum(valid==i_fv)==1)
        if(sum(erode_valid_ids==i_fv)==1)
            final_model(connected_volumes_array(i_fv).IndVoxelList) = 1;
                connected_volumes_array(i_fv).Bone = 1;
        end
    end
    
end


%smoothing step?
binary_map_filled = final_model;
smoothed = imgaussfilt3(double(binary_map_filled),0.4);

smoothed_thresh = smoothed>0.8;
CC_final = bwconncomp(smoothed_thresh);
for i_cc = 1:length(CC_final.PixelIdxList)
    if(length(CC_final.PixelIdxList{i_cc})<500)
        binary_map_filled(CC_final.PixelIdxList{i_cc}) = 0;
    end
end

smoothed = imgaussfilt3(double(binary_map_filled),0.4);

SaveBone3D(size(ct3D),connected_volumes_array, X_ct);

%PLY WRITE
%SaveBone3D(

% comp_map = binary_labelled_map.*2;
% comp_map(final_model==1) = 1;
% figure
% hpatA= ViewSeperateVolumes(final_model,X_ct);
% figure
% hpatB= ViewSeperateVolumes(comp_map,X_ct);
% 
% figure
% col_map = cell2mat(col_arr');
% [hpat] = PATCH_3Darray(fully_labelled_map,X_ct{1},X_ct{2},X_ct{3},col_map,'col');
% 
% 
% for i_va = 1:length(connected_volumes_array)
%     
%     if(connected_volumes_array(i_va).NumVoxels>800)
%        if(isempty(hpat{1,i_va}))
%            i_va
%        else
%         hpat{1,i_va}.FaceAlpha = 0.1;
%         hpat{1,i_va}.EdgeAlpha = 0.1;
%        end
%     else
%         if(isempty(hpat{1,i_va}))
%            i_va
%        else
%         hpat{1,i_va}.FaceAlpha = 0.6;
%         hpat{1,i_va}.EdgeAlpha = 0.6;
%        end
%     end
%     
%   
%     
% end

% 
% yfit = trainedModel.predictFcn(full_table);
% class_predict_map = uint16(zero_map);
% for i_lv = 1:length(valid)
%     if(yfit(i_lv)==0)
%         class_predict_map(connected_volumes(valid(i_lv)).IndVoxelList) = 1;
%     else
%         class_predict_map(connected_volumes(valid(i_lv)).IndVoxelList) = 2;
%     end
% end
% 
% figure
% hpat = ViewSeperateVolumes(class_predict_map, X_ct);



