clear all 
close all
dataset = 1;
%load data
dt = string(Configuration.TorsoDatasets);
dtt = dt(dataset);

[~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(char(dtt), false);


load("C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\PreSeperation\" + num2str(dataset) + ".mat");


[labelled_map] = PreClassifySeperate2(ct3D,binary_map_filled, CT_dimension_spacing,X_ct, false);
    

pre_process_lab_map = labelled_map;
[connected_volumes_array] = PostPreClassifySeparateProcess(pre_process_lab_map, ct3D,dataset,X_ct);

full_label_map = uint8(zeros(size(ct3D)));
for i_cv = 1:length(connected_volumes_array)
    full_label_map(connected_volumes_array(i_cv).IndVoxelList) = i_cv;
end

figure
ViewSeperateVolumes(full_label_map, X_ct);



%load("C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\PostSeperation\" + num2str(dataset) + "_no_coating.mat")

% %post process separate non-connected:
% zero_map = uint8(zeros(size(ct3D)));
% zero_map16 =  uint16(zeros(size(ct3D)));
% 
% 
% initial_num_labels = max(final_labelled_map(:));
% 
% new_label_counter = 1;
% 
% for i_label = 1:initial_num_labels
%     cur_map = zero_map;
%     cur_map(final_labelled_map==i_label) = 1;
%     
%     CC_cur = bwconncomp(cur_map);
%     
%     if(length(CC_cur.PixelIdxList)>1)
%         %adjust labels
%         for i_new_label = 2:length(CC_cur.PixelIdxList)
%             final_labelled_map(CC_cur.PixelIdxList{i_new_label}) = initial_num_labels+new_label_counter;
%             new_label_counter = new_label_counter+1;
%         end
%     end
%     
% end
% 
% new_num_vols = max(final_labelled_map(:));
% %double check:
% d_check = length(unique(final_labelled_map(:))) -1;
% 
%pre_process_lab_map = final_labelled_map;
%[connected_volumes_array, final_labelled_map] = PostPreClassifySeparateProcess(pre_process_lab_map, ct3D,binary_map_filled, CT_dimension_spacing,X_ct);

%PostPreClassifySeparateProcess dev
% full_long_pixelIdx_list = full_long_pixelIdx_list(~cellfun(@isempty, full_long_pixelIdx_list));
% 
% vol_lengths = zeros(new_num_vols,1);
% 
% %of the small vols,
% %   find those without other 'bone' nearby and remove
% for i_label = 1:new_num_vols
%     vol_lengths(i_label) = length(full_long_pixelIdx_list{i_label});
% end
% 
% short_label_min = 100;
% short_labels = find(vol_lengths<short_label_min);
% 
%  rp  = regionprops3(final_labelled_map, {'Centroid','EquivDiameter'});
%  CoM = rp(:,'Centroid');
%  EqDia = rp(:,'EquivDiameter');
% %setup co-ords transform
% %gi = griddedInterpolant(1:size(ct3D,1), 1:size(ct3D,2), 1:size(ct3D,3), ndgrid(X_ct{1},X_ct{2},X_ct{3}));
% gi = {};
% for i = 1:3
%     gi{i} = griddedInterpolant(1:size(ct3D,i), X_ct{i});
% 
% end
% 
% all_indices = cell2mat(full_long_pixelIdx_list);
% 
% 
% small_indices_arr = {};
% 
% for i_short = 1:length(short_labels)
%     cur_label = short_labels(i_short);
%         small_indices_arr{i_short} = full_long_pixelIdx_list{cur_label};
% 
% end
% 
% %remove cur label indices from all indices list
% all_small_indices = cell2mat(small_indices_arr');
% other_indices = setdiff(all_indices, all_small_indices);
% 
% %use NN tagging to find distance to closest part of the model
% 
% %unlabelled_point_indices = setdiff(unlabelled_point_indices,all_seed_indexes);
% 
% [ul_arr] = CoM{short_labels, 'Centroid'};
% 
% %labelled_point_indices = find(labelled_map);
% %labelled_point_indices = all_seed_indexes;
% [lp_1,lp_2,lp_3] = ind2sub(size(final_labelled_map),other_indices);
% %convert to real positions using x_ct
% ulp_1 = gi{1}(ul_arr(:,2));%swapped on purpose!
% ulp_2 = gi{2}(ul_arr(:,1));%
% ulp_3 = gi{3}(ul_arr(:,3));
% 
% lp_1 = X_ct{1}(lp_1);
% lp_2 = X_ct{2}(lp_2);
% lp_3 = X_ct{3}(lp_3);
% 
% 
% 
% [~, D] = knnsearch([lp_1',lp_2',lp_3'],[ulp_1,ulp_2,ulp_3]);
% 
% 
% f = figure
% binary_labelled_map = zero_map;
% binary_labelled_map(all_indices) = 1;
% 
% short_rs =  EqDia{short_labels, 'EquivDiameter'};
% short_rs = short_rs/2;
% 
% D_min_r = D-short_rs;
% 
% %dist_range = 1:80;
% min_allowed_distance = 3;
% dist_inds = find(D_min_r>=min_allowed_distance); %to_purge
% %dist_inds = find(D_min_r<3); % survigin
% culled_labels = short_labels(dist_inds);
% % target_vol = sel_labels;
% % use_indices = true;
% % target_vol_inds = {};
% % for i = 1:length(sel_labels)
% %     target_vol_inds{i} = full_long_pixelIdx_list{sel_labels(i)};
% % end
% 
% 
% %view_individual_volume_script;
% valid_volumes = setdiff(1:length(full_long_pixelIdx_list),culled_labels);
% 
% 
% %split volumes into large, medium and small for easier classification!
% valid_lengths = vol_lengths(valid_volumes);
% 
% v_small_boundary = 10;
% v_small_vols = valid_volumes(find(valid_lengths<=v_small_boundary));
% 
% small_boundary = 50;
% small_vols = valid_volumes(find(valid_lengths<=small_boundary));
% only_small_vols = setdiff(small_vols, v_small_vols);
% 
% 
% 
% 
% medium_boundary = 200;
% medium_vols = valid_volumes(find(valid_lengths<=medium_boundary));
% only_medium_vols = setdiff(medium_vols, small_vols);
% 
% large_boundary = 500;
% large_vols = valid_volumes(find(valid_lengths<=large_boundary));
% only_large_vols = setdiff(large_vols, medium_vols);
% 
% 
% v_large_vols = valid_volumes(find(valid_lengths>large_boundary));
% only_v_large_vols = setdiff(v_large_vols, large_vols);
% 
% %make a label map by vol size for good picture
% size_label_map = zero_map;
% 
% 
% cur_list = v_small_vols;
% cur_lab = 1;
% for i_v = 1:length(cur_list)
%     size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
% end
% 
% cur_list = only_small_vols;
% cur_lab = 2;
% for i_v = 1:length(cur_list)
%     size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
% end
% 
% cur_list = only_medium_vols;
% cur_lab = 3;
% for i_v = 1:length(cur_list)
%     size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
% end
% 
% cur_list = only_large_vols;
% cur_lab = 4;
% for i_v = 1:length(cur_list)
%     size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
% end
% 
% cur_list = only_v_large_vols;
% cur_lab = 5;
% for i_v = 1:length(cur_list)
%     size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
% end
% 
% cur_list = culled_labels;
% cur_lab = 6;
% for i_v = 1:length(cur_list)
%     size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
% end
% 
% f = figure
% [hpat] = ViewSeperateVolumes(size_label_map, X_ct);
% 
% for i_h = 1:6
%     hpat{i_h}.FaceAlpha = 0.1;
%     hpat{i_h}.EdgeAlpha = 0.1;
% end

%save("C:\Users\mazna\Documents\nl\U\P\Report\FinalReport\FigureData\HardSeperate\final_" + num2str(dataset),'f');

% %view data
% figure
% ViewSeperateVolumes(final_labelled_map, X_ct);
% %
% 
% figure
% imshow3D(final_labelled_map);
% 
% 
% %through spine view:
% x1= 253; y1 = 298;
% z1 = 1:135;
% x2 = 255; y2 = 259;
% z2 = 136:200;
% 
% spine_part_1 = squeeze(final_labelled_map(y1,x1,z1));
% 
% spine_part_2 = squeeze(final_labelled_map(y2,x2,z2));
% figure
% hold on
% plot(squeeze(spine_part_1))
% plot(squeeze(spine_part_2))
% 
% f = figure;
% 
% 
% binary_labelled_map = zero_map;
% binary_labelled_map(final_labelled_map~=0) = 1;
% 
% %view specific vol
% 
% target_vol = [106, 105, 100, 114, 113, 115];
% 
% view_individual_volume_script;

% 
% spec_map = zero_map;
% 
% spec_map(final_labelled_map==target_vol) = 1;
% 
% %change bounds to focus on target vol
% 
% target_indices = find(final_labelled_map==target_vol);
% [target_sub1, target_sub2, target_sub3] = ind2sub(size(final_labelled_map),target_indices);
% 
% border_radius = 40;
% 
% min_x1 = max(1, min(target_sub1)-border_radius);
% min_x2 = max(1, min(target_sub2)-border_radius);
% min_x3 = max(1, min(target_sub3)-border_radius);
% 
% max_x1 = min(size(ct3D,1), max(target_sub1)+border_radius);
% max_x2 = min(size(ct3D,2), max(target_sub2)+border_radius);
% max_x3 = min(size(ct3D,3), max(target_sub3)+border_radius);
% 
% 
% spec_map = spec_map(min_x1:max_x1, min_x2:max_x2, min_x3:max_x3);
% 
% pos = get(f, 'Position'); %// gives x left, y bottom, width, height
% close(f);
% 
% f = figure('Name', "cur_target_vol: " + num2str(target_vol));
% set(f, 'Position', pos);
% figure(f)
% 
% %clf(f)
% hold  on
% [hpat1] = ViewSeperateVolumes(binary_labelled_map(min_x1:max_x1, min_x2:max_x2, min_x3:max_x3), X_ct, min_x1:max_x1, min_x2:max_x2, min_x3:max_x3);
% hpat1.FaceAlpha  = 0;
% hpat1.EdgeAlpha  = 0.1;
% [hpat2] = ViewSeperateVolumes(spec_map, X_ct, min_x1:max_x1, min_x2:max_x2, min_x3:max_x3);
% hpat2.FaceColor = 'g';
% hpat2.FaceAlpha  = 0.2;
% hpat2.EdgeColor = 'g';
% hpat2.EdgeAlpha  = 0.2;
% save("C:\Users\mazna\Documents\nl\U\P\Report\FinalReport\FigureData\HardSeperate\IndVolumeAnalysis\" + num2str(target_vol), 'f');