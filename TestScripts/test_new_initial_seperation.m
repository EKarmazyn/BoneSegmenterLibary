%first test DB 1
clear all
close all
dbstop if error
dataset = 1;
tic
dt = string(Configuration.TorsoDatasets);
dtt = dt(dataset);
[initial_model, X_ct, minimum_bone_intensity] = SegmentBoneV2(char(dtt), 'end_after_fuzzy_find_bone', 1);
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
%remove bed (by sampling 20 points from each, and checking flatness of
%surroundings
% test_r = 10;
% invalid_r = test_r+5;
% szModel = size(labelled_initial_model);
% bed_remove = [];
% for i_CC = indices_CC
%     points = find(labelled_initial_model==i_CC);
%     [stemp{1},stemp{2},stemp{3}] = ind2sub(size(labelled_initial_model),points);
%     invalid = (stemp{1}<invalid_r) + (stemp{1}>(szModel(1)-invalid_r)) ...
%                 + (stemp{2}<invalid_r) + (stemp{2}>(szModel(3)-invalid_r)) ...
%                 + (stemp{3}<invalid_r) + (stemp{3}>(szModel(3)-invalid_r)) ;
%     
%     invalid(invalid>0) = 1;
%     points(invalid==1) = [];
%     if(isempty(points))
%         continue;
%     end
%     sample_points = round((0.01:0.05:0.99)*length(points));
%     results_store = [];
%     for i_sp = 1:length(sample_points)
%         %extract region
%         [s1,s2,s3] = ind2sub(szModel,points(sample_points(i_sp)));
%         sample_region = labelled_initial_model(s1-test_r:s1+test_r,s2-test_r:s2+test_r, s3-test_r:s3+test_r);
%         %crop out planes
%         
%         
%         %remove zero
%             plane1_toRemove = [];
%             for iter_plane1 = 1:size(sample_region,3)
%                 if(~any(sample_region(:,:,iter_plane1)))
%                     plane1_toRemove = [plane1_toRemove,  iter_plane1];
%                 end
%             end
%             szRes = size(sample_region);
%             for i_sz = 3:-1:(length(szRes)+1)
%                 szRes(i_sz) = 1;
%             end
%             szRes(3) = szRes(3)-length(plane1_toRemove);
%             sample_region(:,:,plane1_toRemove) = [];
%             sample_region = reshape(sample_region,szRes);
%             
%             plane2_toRemove = [];
%             for iter_plane2 = 1:size(sample_region,1)
%                 if(~any(sample_region(iter_plane2,:,:)))
%                     plane2_toRemove = [plane2_toRemove,  iter_plane2];
%                 end
%             end
%             szRes = size(sample_region);
%             for i_sz = 3:-1:(length(szRes)+1)
%                 szRes(i_sz) = 1;
%             end
%             szRes(1) = szRes(1)-length(plane2_toRemove);
%             sample_region(plane2_toRemove,:,:) = [];
%             sample_region = reshape(sample_region,szRes);
%             
%             plane3_toRemove = [];
%             for iter_plane3 = 1:size(sample_region,2)
%                 if(~any(sample_region(:,iter_plane3,:)))
%                     plane3_toRemove = [plane3_toRemove,  iter_plane3];
%                 end
%             end
%             szRes = size(sample_region);
%             for i_sz = 3:-1:(length(szRes)+1)
%                 szRes(i_sz) = 1;
%             end
%             szRes(2) = szRes(2)-length(plane3_toRemove);
%             sample_region(:,plane3_toRemove,:) = [];
%             sample_region = reshape(sample_region,szRes);
%             
%             
%             szRes = size(sample_region);
%             for i_sz = 3:-1:(length(szRes)+1)
%                 szRes(i_sz) = 1;
%             end
%             
%             results_store = [results_store; szRes];
%             
%     end
%     avgs = mean(results_store);
%     if(any(avgs<5))
%         bed_remove = [bed_remove; i_CC];
%     end
% end
% 
% for i_bedRemove = 1:length(bed_remove)
%     labelled_initial_model(CC.PixelIdxList{i_bedRemove}) = 0;
% end

%generate edge map
dog_gf1 = 0.0104;
dog_gf2 = 1.9997;
dog_filt = imgaussfilt3(double(ct3D),dog_gf1)-imgaussfilt3(double(ct3D),dog_gf2);
tic
smoothed = imgaussfilt3(dog_filt,4,'FilterSize',11);
toc

dog_minus_smoothed = dog_filt - smoothed;

basic_surfaces = dog_minus_smoothed>50;
CC_surf = bwconncomp(basic_surfaces);
for i_CC = 1:length(CC_surf.PixelIdxList)
    if(length(CC_surf.PixelIdxList{i_CC}) < 100)
        basic_surfaces(CC_surf.PixelIdxList{i_CC}) = 0;
        
    end
end

% figure
% lengths = [];
% for i_CC = 1:length(CC_surf.PixelIdxList)
%     lengths(i_CC) = length(CC_surf.PixelIdxList{i_CC});
% end
% histogram(lengths, 1:length(CC_surf.PixelIdxList));


%basic_surfaces(initial_model==0) = 0;

%for each point in the 'basic_surfaces', calculate nearby regionprops
tic
[Gx,Gy,Gz] = imgradientxyz(smoothed);
toc


%dil_r = 4;
dil_r = 2;

Gx(basic_surfaces==0) = 0;
Gy(basic_surfaces==0) = 0;
Gz(basic_surfaces==0) = 0;
Gx2 = Gx.^2;
Gy2 = Gy.^2;
Gz2 = Gz.^2;
Gsum = Gx2+Gy2+Gz2;
g_x_ratios = round((Gx2./Gsum)*dil_r);
g_y_ratios = round((Gy2./Gsum)*dil_r);
g_z_ratios = round((Gz2./Gsum)*dil_r);
g_x_ratios(basic_surfaces==0) = 0;
g_y_ratios(basic_surfaces==0) = 0;
g_z_ratios(basic_surfaces==0) = 0;
gRatiosCell = {g_x_ratios, g_y_ratios, g_z_ratios};

%filter shapes (start with sphere)

sphere_SE = strel('sphere',dil_r);
sphere_patch = sphere_SE.Neighborhood;
%x/y/z_stretch = 0:4 where 4 gives full sphere width and 0 gives only the
%central pixels
%we have ratios which sum up to 1;

%we have 5^3 different filters
filter_bank = cell(dil_r*3,1);
filter_bank_dir_values = cell(dil_r*3,1);
filter_iterator = 1;
for i_x = 1:dil_r+1
    for j_y = 1:dil_r+1
        for k_z = 1:dil_r+1
            if(i_x + j_y + k_z ~= dil_r+3)
                continue;
            end
            %iter = (i_x-1)*5^2 + (j_y-1)*5 + k_z;
            szRes = [(i_x-1)*2 + 1, (j_y-1)*2 + 1, (k_z-1)*2 + 1];
            resized = imresize3(uint8(sphere_patch), szRes);
            resized = reshape(resized,szRes);
            %remove zero
            plane1_toRemove = [];
            for iter_plane1 = 1:size(resized,3)
                if(~any(resized(:,:,iter_plane1)))
                    plane1_toRemove = [plane1_toRemove,  iter_plane1];
                end
            end
            szRes = size(resized);
            for i_sz = 3:-1:(length(szRes)+1)
                szRes(i_sz) = 1;
            end
            szRes(3) = szRes(3)-length(plane1_toRemove);
            resized(:,:,plane1_toRemove) = [];
            resized = reshape(resized,szRes);
            
            plane2_toRemove = [];
            for iter_plane2 = 1:size(resized,1)
                if(~any(resized(iter_plane2,:,:)))
                    plane2_toRemove = [plane2_toRemove,  iter_plane2];
                end
            end
            szRes = size(resized);
            for i_sz = 3:-1:(length(szRes)+1)
                szRes(i_sz) = 1;
            end
            szRes(1) = szRes(1)-length(plane2_toRemove);
            resized(plane2_toRemove,:,:) = [];
            resized = reshape(resized,szRes);
            
            plane3_toRemove = [];
            for iter_plane3 = 1:size(resized,2)
                if(~any(resized(:,iter_plane3,:)))
                    plane3_toRemove = [plane3_toRemove,  iter_plane3];
                end
            end
            szRes = size(resized);
            for i_sz = 3:-1:(length(szRes)+1)
                szRes(i_sz) = 1;
            end
            szRes(2) = szRes(2)-length(plane3_toRemove);
            resized(:,plane3_toRemove,:) = [];
            resized = reshape(resized,szRes);
            
            %make square
            max_dim = max(size(resized));
            square = zeros(max_dim,max_dim,max_dim);
            szRes = size(resized);
            for i_sz = 3:-1:(length(szRes)+1)
                szRes(i_sz) = 1;
            end
            offset = max_dim*ones(3,1)' - szRes;
            offset = offset./2;
            
            
           
            square(1+offset(1):end-offset(1),1+offset(2):end-offset(2),1+offset(3):end-offset(3)) = resized;
            
            filter_bank{filter_iterator} = square;
            filter_bank_dir_values{filter_iterator} = [i_x-1, j_y-1, k_z-1];
            filter_iterator = filter_iterator+1;
        end
    end
end

g_filter_i = zeros(size(g_x_ratios));
gCell = cell(3,dil_r+1);
for i_dim = 1:3
    for i_gRatio = 1:dil_r+1
        gCell{i_dim,i_gRatio} = (gRatiosCell{i_dim} == (i_gRatio-1));
    end
end

for i_filter = 1:length(filter_bank_dir_values)
    g_filter_i((gCell{1,filter_bank_dir_values{i_filter}(1)+1} .* gCell{2,filter_bank_dir_values{i_filter}(2)+1} .* gCell{3,filter_bank_dir_values{i_filter}(3)+1}) == 1) = i_filter;
    
end
g_filter_i(basic_surfaces==0) = 0;





dilated_surfaces_map = basic_surfaces;
for i_filter = 1:length(filter_bank_dir_values)
    tic
    cur_surf_points = (g_filter_i==i_filter);
    dil_cur_surf_points = imdilate(cur_surf_points, filter_bank{i_filter});
    dilated_surfaces_map(dil_cur_surf_points==1) = 1;
    toc
end

near_bone_mask =  imdilate(initial_model, [1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1]);
dilated_surfaces_map_near_bone = dilated_surfaces_map;
dilated_surfaces_map_near_bone(near_bone_mask==1) = 1;


%find internal bone
summed_internal_bone = zeros(size(initial_model));
input_mask = near_bone_mask;
initial_bone_model_plus_additions = initial_model;
for i_internalPoits = 1:5
    [initial_bone_model_plus_additions,internal_map, num_dir_bone_found_in] = InternalPointsInline(input_mask,initial_bone_model_plus_additions,30,1);
    summed_internal_bone(internal_map==1) = 1;
    
    input_mask = imdilate(internal_map, [1 1 1 ; 1 1 1 ; 1 1 1 ]);
    input_mask(initial_bone_model_plus_additions==1) = 0;
end
input_mask = imdilate(summed_internal_bone, [1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1]);
input_mask(summed_internal_bone==1) = 0;
[updated_bone_map,internal_map, num_dir_bone_found_in] = InternalPointsInline(input_mask,initial_bone_model_plus_additions,10,2);
summed_internal_bone(internal_map==1) = 1;
initial_bone_model_plus_additions = updated_bone_map;


clear Gx Gx2 Gy Gy2 Gz Gz2 Gsum gCell gRatiosCell

%assign summed_internal_bone to ther labelled_dilated_model




%labelled_dilated_model = labelmatrix(CC);
zero_map = zeros(size(labelled_dilated_model));
new_labelled_initial_model = uint16(zero_map);
%new_model = new_labelled_initial_model;
%new_model(labelled_dilated_model>0) = 1;
new_model(summed_internal_bone==1) = 1;
new_model_minus_surfaces = new_model;
new_model_minus_surfaces(dilated_surfaces_map==1) = 0;
new_label_counter = 1;
new_model_cc = bwconncomp(new_model_minus_surfaces);
labelled_new_model = labelmatrix(new_model_cc);
for initial_cv = indices_CC
    %seperate the dilated CV areas then CC analysis and assign labels to
    %initial model 
    cur_cv = zero_map;
    cur_cv(CC.PixelIdxList{initial_cv}) = 1;
    cur_cv(dilated_surfaces_map==1) = 0;
    CC_specific = bwconncomp(cur_cv);
    new_label_matrix = uint16(labelmatrix(CC_specific));
    new_labelled_initial_model = new_labelled_initial_model+(new_label_matrix+(new_label_counter-1));
    new_label_counter = new_label_counter+length(CC_specific.PixelIdxList);
end
    



