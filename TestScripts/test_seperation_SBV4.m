


xy_thresh = 15
z_thresh = 0
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
%figure
%imshow3D(updated_binary_map);

%seperate any non-connected volumes!
% 
%such as slow way of checking... it seems always fine
% for i_label = 1:length(labels)
%     tic
%     CC = bwconncomp(updated_binary_map(fully_labelled==labels(i_label)));
%     toc
%     if(length(CC.PixelIdxList) > 1)
%         display('error');
%     end
% end

%get index list for each label

labelled_indices = find(updated_binary_map);
indices_labels = fully_labelled(labelled_indices);

label_index_list = cell(length(labels),1);
tic
for i_label = 1:length(labels)
     %tic
        label_index_list{i_label,1} = labelled_indices(indices_labels==labels(i_label));
     %toc
end
toc

valid_labelled_index_list = cell(length(label_index_list)-1,1);
tic
for i_label = 1:length(labels)-1
     %tic
        valid_labelled_index_list{i_label,1} = label_index_list{i_label+1};
     %toc
end
toc



tic
connected_volumes_array = ConnectedVolume.FromLabelledIndices(valid_labelled_index_list);
toc
%pick out some error volumes and observe properties

%mark connected volumes with bone seeds
tic
[prob_bone_seed_map, prob_minimum_bone_seed_intensity, def_bone_seed_map, def_bone_seed_min_intensity] = FindBoneSeeds(ct3D, ones(size(ct3D)));


ConnectedVolume.MarkWithBoneSeeds(connected_volumes_array,prob_bone_seed_map, def_bone_seed_map);
t = toc;
fprintf("MarkWithBoneSeeds  Time-Taken: " + num2str(t)+"\n");

%extract features of connected volumes
T = tic;
for cv = 1:length(connected_volumes_array)
    %t_single_start = tic;
    connected_volumes_array(cv).GenerateIndividualFeatures(ct3D);
    %t = toc(t_single_start);
    %fprintf("single_cv_features  Time-Taken: " + num2str(t)+"\n");
end

t = toc(T);
fprintf("all_cv_features  Time-Taken: " + num2str(t)+"\n");




T = tic;
ConnectedVolume.GenerateRegionPropsFeatures(connected_volumes_array, fully_labelled);
t = toc(T);
fprintf("GenerateRegionPropsFeatures  Time-Taken: " + num2str(t)+"\n");

% twenty_only = zeros(size(updated_binary_map));
% twenty_only(fully_labelled==20) = 1;

%%%%%%%%%%%%%%%%%%% new strat:

only_large_labelled_map = fully_labelled;
for i_label = 1:length(valid_labelled_index_list)
     if(length(valid_labelled_index_list{i_label,1})<100)
         only_large_labelled_map(valid_labelled_index_list{i_label,1}) = 0;
     end
end

figure
ViewSeperateVolumes(only_large_labelled_map,X_ct);
figure
imshow3D(only_large_labelled_map);

%%%%%%%%%%%%%%%%%%%%555

error_labels = [20,677,679,5918,5920,5922];
for i_error = 1:length(error_labels)
    connected_volumes_array(error_labels(i_error)).Bone = -1;
    
end

bone_labels = [81, 14, 7259, 7308, 1367];
for i_bone = 1:length(bone_labels)
    connected_volumes_array(bone_labels(i_bone)).Bone = 1;
    
end

T = tic;
[full_table, bone_table, non_bone_table, labelled_table, unlabelled_table] = ConnectedVolume.ConstructDataTable(connected_volumes_array);
t = toc(T);
fprintf("ConstructDataTable  Time-Taken: " + num2str(t)+"\n");

rows = full_table.BoneType==0;
extended_table = full_table;
extended_table.BoneType(rows) = 1;

for i = 1:100
    extended_table = [extended_table; non_bone_table];
end


%predict
%yfit = fineGaussSVMModel.predictFcn(full_table);
yfit = coarseTreeModel.predictFcn(full_table);


bone_ind = (yfit==1);
not_bone_ind = (yfit==-1);

zero_map = zeros(size(ct3D));

bone_vols = connected_volumes_array(bone_ind);

bone_map = zero_map;


points_ind_list = [];
for i_bi = 1:length(bone_vols)
    cur_co = bone_vols(i_bi);
    
    points_ind_list = [points_ind_list; cur_co.IndVoxelList];
    
end

bone_map(points_ind_list) = 1;


not_bone_vols = connected_volumes_array(not_bone_ind);
not_bone_map = zero_map;

points_ind_list = [];
for i_bi = 1:length(not_bone_vols)
    cur_co = not_bone_vols(i_bi);
    
    points_ind_list = [points_ind_list; cur_co.IndVoxelList];
    
end

not_bone_map(points_ind_list) = 1;


figure('Name', 'bone')
imshow3D(bone_map,[0 1])

figure('Name', 'not_bone_map')
imshow3D(not_bone_map,[0 1])
























% 
% 
% dbstop if error
% 
% for dataset = 1:10
%     
%     
%     dt = string(Configuration.TorsoDatasets);
%     dtt = dt(dataset);
%     
%     
%     %load
%     load("C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\PreSeperation\" + num2str(dataset) + ".mat");
%     [~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(char(dtt), false);
% 
%     %report
%     start_time = tic;
%     
%     % --> binary_map_filled <--
%     
%     %voided = HardSeperateVolumes(ct3D,binary_map_filled, CT_dimension_spacing, 0,10);
%     for xy_thresh = -5:5:5
%         for z_thresh = 0:5:10
%             voided = HardSeperateVolumes(ct3D,binary_map_filled, CT_dimension_spacing, xy_thresh,z_thresh);
%             [connected_volumes2, labelled_binary_map2,all_points_ind_list2] = ConnectedVolumes(voided);
%             figure
%             xy_thresh
%             z_thresh
%             imshow3D(labelled_binary_map2);
%         end
%     end
%     
% %     %test SeperateVolumesV6
% %     
% %     [connected_volumes, labelled_binary_map,all_points_ind_list] = ConnectedVolumes(binary_map_filled);
% %     tElapsed = toc(start_time);
% %     fprintf("Time-Taken-CV-Generation: " + num2str(tElapsed)+"\n");
% %     
% %     start_time2 = tic;
% %     
% %     dog_gf1 = 0.0104;
% %     dog_gf2 = 1.9997;
% %     dog_filt_threshold = 0;
% %     seperation_use_model_surface = false;
% %     seperation_min_voxels = 6;
% %     dilate_join_internals = true;
% %     
% %     
% %     dog_filt = imgaussfilt3(double(ct3D),dog_gf1)-imgaussfilt3(double(ct3D),dog_gf2);
% % 
% %     
% %     smoothed = imgaussfilt3(dog_filt,4,'FilterSize',11);
% %     
% %     
% %     dog_minus_smoothed = dog_filt - smoothed;
% %     
% %     
% %     
% %     
% %     
% %     bt = 0.3;
% %     span = 2;
% %     sps = 4;
% %     h = gaussdesign(bt,span,sps);
% %     %fvtool(h,'impulse')
% %     
% %     
% %     span = 4;
% %     sps = 2;
% %     h2 = gaussdesign(bt,span,sps);
% %     %fvtool(h2,'impulse')
% %     
% %     h_dog = h2-h;
% %     %fvtool(h_dog,'impulse')
% %     
% %     %h_dog2 = h_dog(11:23);
% %     h_dog2 = h_dog;
% %     
% %     filt_length_1d = length(h_dog1);
% %     
% %     cur_test_filt = h_dog2;
% %     
% %     
% %     
% %     %edge X (ACTUALLY UPDOWN (i,e. along cols in imshow3d)
% %     x_edge_map = zeros(size(ct3D));
% %     offset_raw = int32(ct3D) - 1000;
% %     
% %     for i_y = 1:512
% %         %tic
% %         for i_z = 1:361
% %             x_array = offset_raw(:,i_y,i_z);
% %             edge_v = conv(x_array,cur_test_filt,'same');
% %             x_edge_map(:,i_y,i_z) = edge_v(:);
% %         end
% %         %toc
% %     end
% %     
% %     figure('Name','x')
% %     imshow3D(x_edge_map, [-500,500]);
% %     
% %     
% %     %edge Y
% %     
% %     
% %     y_edge_map = zeros(size(ct3D));
% %     %offset_raw = int32(ct3D) - 1000;
% %     
% %     for i_x = 1:512
% %         %tic
% %         for i_z = 1:361
% %             x_array = squeeze(offset_raw(i_x,:,i_z));
% %             edge_v = conv(x_array,cur_test_filt,'same');
% %             y_edge_map(i_x,:,i_z) = edge_v(:);
% %         end
% %         %toc
% %     end
% %     figure('Name','y')
% %     imshow3D(y_edge_map, [-500,500]);
% %     
% %     
% %     %xy_sum_edge_map = zeros(size(ct3D));
% %     xy_sum_edge_map = (x_edge_map.^2 + y_edge_map.^2).^(1/2);
% %     
% %     
% %     %reverse DoG filter to find double edges? - idea but maybe bad is
% %     %'hard' areas
% %     
% %     
% %     
% %     
% %     
% %     
% %     %figure('Name','xy')
% %     %imshow3D(xy_sum_edge_map, [-500,500]);
% %     
% %     dog_gf1 = 0.0104;
% %     dog_gf2 = 1.9997;
% %     xy_dog_sum_edge_map = zeros(size(ct3D));
% %     for i_z = 1:361
% %         xy_dog_sum_edge_map(:,:,i_z) = imgaussfilt(xy_sum_edge_map(:,:,i_z),dog_gf1) - ...
% %             imgaussfilt(xy_sum_edge_map(:,:,i_z),dog_gf2);
% %      end
% %     
% %     figure('Name','xy_dog_sum_edge_map')
% %     imshow3D(xy_dog_sum_edge_map, [0,1]);
% %     xy_dog_sum_edge_map_culled = xy_dog_sum_edge_map;
% %     xy_dog_sum_edge_map_culled(ct3D<1000) = -100;
% %     
% %     
% %     binary_xy_edge = xy_dog_sum_edge_map_culled>0;
% %     figure('Name','binary_xy_edge')
% %     imshow3D(binary_xy_edge, [0,1]);
% %     
% %     
% %     
% %     
% %     %%% TESTING dog of edge vs edge
% %     figure('Name','binary_map_filled')
% %     imshow3D(binary_map_filled, [0,1]);
% %     
% %     xy_binary_map = zeros(size(xy_sum_edge_map));
% %     xy_binary_map(xy_sum_edge_map>50) = 1;
% %     xy_binary_map20 = zeros(size(xy_sum_edge_map));
% %     xy_binary_map20(xy_sum_edge_map>20) = 1;
% %     
% %     dog_xy_binary_map = zeros(size(xy_sum_edge_map));
% %     dog_xy_binary_map(xy_dog_sum_edge_map>0) = 1;
% %     
% %     
% %     cur_test_void = dog_xy_binary_map;
% %     voided = binary_map_filled;
% %     voided(cur_test_void==1) = 0;
% %     %figure('Name','dog_xy_binary_map0_voided')
% %     %imshow3D(voided, [0 1]);
% %     
% %     %xy_binary_map20 
% %     
% %     for i_z = 1:361
% %         voided(:,:,i_z) = imerode(voided(:,:,i_z), strel('diamond',1));
% %     end
% %     
% %     for i_z = 1:361
% %         voided2(:,:,i_z) = imdilate(voided(:,:,i_z), strel('diamond',1));
% %     end
% %     %figure('Name','diltest')
% %     %imshow3D(voided2, [0 1]);
% %     
% %     
% %     
% %     %edge Z %%%%%%%%%%%%%%%%%%%%%%%%%%5
% %     %dog_xy_binary_map ,__ already got xy
% %     
% %     %adjust filter
% %     z_rel_scale = CT_dimension_spacing(3)/CT_dimension_spacing(1);
% %     
% %     bt = 0.3;
% %     span = 1;
% %     sps = 2;
% %     hz = gaussdesign(bt,span,sps);
% %     
% %     span = 2;
% %     sps = 1;
% %     h2z = gaussdesign(bt,span,sps);
% %     hz_dog = h2z-hz;
% %     fvtool(hz_dog,'impulse')
% %     
% %     
% %     z_edge_map = zeros(size(ct3D));
% %     for i_x = 1:512
% %         %tic
% %         for i_y = 1:512
% %             z_array = squeeze(offset_raw(i_x,i_y,:));
% %             edge_v = conv(z_array,cur_test_filt,'same');
% %             z_edge_map(i_x,i_y,:) = edge_v(:);
% %         end
% %         %toc
% %     end
% %     
% %     %xyz_sum_edge_map = (x_edge_map.^2 + y_edge_map.^2+ z_edge_map.^2).^(1/2);
% %     
% %     %figure('Name','z_edge_map')
% %     %imshow3D(z_edge_map, [0,1]);
% %     
% %     z_binary_map = zeros(size(xy_sum_edge_map));
% %     z_binary_map(z_edge_map>10) = 1;
% %   
% %     
% %     
% %     voided = binary_map_filled;
% %     voided(dog_xy_binary_map==1) = 0;
% %     voided(z_binary_map==1) = 0;
% %      
% %     figure('Name','voided')
% %     imshow3D(voided, [0,1]);
% %     
% %     start_time = tic;
% %     
% %     [connected_volumes2, labelled_binary_map2,all_points_ind_list2] = ConnectedVolumes(voided);
% %     tElapsed = toc(start_time);
% %     fprintf("Time-Taken-CV-Generation2: " + num2str(tElapsed)+"\n");
% %     
% %     %figure
% %     %ViewSeperateVolumes(labelled_binary_map2, X_ct);
% %     
% %     
%   
% %erode
% 
% %dilate
% 
% %cc
% 
% %coating
% 
% %knn labelling
%     
% %     dog_gf1 = 0.0104;
% %     dog_gf2 = 1.9997;
% %     xy_dog_sum_edge_map = zeros(size(ct3D));
% %     for i_z = 1:361
% %         xy_dog_sum_edge_map(:,:,i_z) = imgaussfilt(xy_sum_edge_map(:,:,i_z),dog_gf1) - ...
% %             imgaussfilt(xy_sum_edge_map(:,:,i_z),dog_gf2);
% %      end
% %     
% %     
% %     xy_dog_sum_edge_map_culled = xy_dog_sum_edge_map;
% %     xy_dog_sum_edge_map_culled(ct3D<1000) = -100;
% %     
% %     
% %     binary_xy_edge = xy_dog_sum_edge_map_culled>0;
%     
%     
%     
%     
%     
%     
%     
%     
%     
% %     tic
% %     %canny_edge_map = canny(double(squeeze(ct3D(:,:,1))));
% %     canny_edge_map = canny(double(ct3D)); 
% %     toc
%     
% %     %edge X
% %     x_dog_filt = reshape(cur_test_filt, filt_length_1d, 1,1 );
% %     x_edges = imfilter(int32(ct3D),x_dog_filt);
% %     figure('Name','x')
% %     imshow3D(x_edges, [-500,500]);
% %     %edge Y
% %     y_dog_filt = reshape(cur_test_filt, 1,filt_length_1d, 1 );
% %     %y = imfilter(int32(ct3D),y_dog_filt,'conv');
% %     y = convn(int32(ct3D),y_dog_filt,'same');
% %     figure('Name','y')
% %     imshow3D(y_edges);
% %     %edge Z
% %     z_dog_filt = reshape(cur_test_filt, 1,1,filt_length_1d );
% %     z_edges = imfilter(int32(ct3D),z_dog_filt,);
% %     figure('Name','z')
% %     imshow3D(z_edges);
%     
%     %surfaces:
%     %either ct3d>1250   or  ct3d>1100 & slight edge(either across slice or
%     %through slice)
%     
%     
%     
%     
%     
%     %report
%     tElapsed = toc(start_time2);
%     fprintf("Time-Taken-Seperate: " + num2str(tElapsed)+"\n");
%     
%     %view
%     figure;
%     imshow3D(new_labelled_map, [0, max(new_labelled_map(:))]);
%     figure;
%     ViewSeperateVolumes(new_labelled_map, X_ct)
%     
%     %hopefully dispose 
%     clear all
%     close all
% end