% 
% dataset = 1;
% tic
% dt = string(Configuration.TorsoDatasets);
% dtt = dt(dataset);
% [new_labelled_map, X_ct, minimum_bone_intensity] = SegmentBoneV2(char(dtt));
% toc
% % 
% 
% %observe distribution of CC
figure
lengths = [];
for i_CC = 1:length(CC.PixelIdxList)
    lengths(i_CC) = length(CC.PixelIdxList{i_CC});
end
histogram(lengths, 1:length(CC.PixelIdxList));
% 
% %observe cutoff at difference lengths of CC
% 
% specific_cutoff = cur_model;
% rem_points = cur_model;
% rem_points = uint8(rem_points);
% CC = bwconncomp(specific_cutoff);
% for i_cc = 1:length(CC.PixelIdxList)
%     if(length(CC.PixelIdxList{i_cc})<100)
%         specific_cutoff(CC.PixelIdxList{i_cc}) = 0;
%         rem_points(CC.PixelIdxList{i_cc}) = 3;
%     end
% end
% image = rem_points;
% figure
% imshow3D(image,[0 max(max(max(image)))]);
% 
% 
clear all
%load workspace
% %load('G:\COverflowDump\nl_\U\P\Data\PostPerforationFixWorkspaces\orig_perf_regions.mat');
 load('G:\COverflowDump\nl_\U\P\Data\PostPerforationFixWorkspaces\final_.mat');

%  
%  %fix_regions_mask_wo_bone
%  
% %  image = fix_regions_mask_wo_bone;
% %  image(cur_model==1) = 5;
% % figure
% % imshow3D(image,[0 max(max(max(image)))]);
% %  
%  %test_mask = zeros(size(cur_model));
%  test_mask = imdilate(cur_model,[1,1,1;1,1,1;1,1,1]);
%  test_mask(cur_model==1) = 0;
dbstop if error

test_mask = imdilate(cur_model,[1,1,1;1,1,1;1,1,1]);
test_summation = test_mask;
test_mask(cur_model==1) = 0;
summed_bone_map = cur_model;
figure
hold on;
for i = 1:20
    [updated_bone_map,internal_map, num_dir_bone_found_in{i}] = InternalPointsInline(test_mask,cur_model,10, 2);
   test_mask = imdilate(internal_map,[1,1,1;1,1,1;1,1,1]);
   test_mask(test_summation==1) = 0;
   test_summation(test_mask==1) = 1;
   summed_bone_map(internal_map==1) = 1;
   
   image = summed_bone_map;
   image(test_mask==1)=3;
   image(cur_model==1)=2;
   
   image(internal_map == 1) = 4;
   subplot(4,5,i)
   imshow(image(209:310,95:160,295),jet(5));
   title(num2str(i));
   
end

%vis code to see which dir is wrong
figure
hold on;
for i = 1:14
    subplot(4,4,i)
    image = squeeze(nearby_bone_map_by_direction(i,:,:,:));
    image(input_mask==0) = 2;
    imshow(image(209:310,95:160,295),jet(3));
end




figure
image = num_dir_bone_found_in{5};
image(cur_model==1)=100;
imshow3D(image,[0 max(max(max(image)))]);



