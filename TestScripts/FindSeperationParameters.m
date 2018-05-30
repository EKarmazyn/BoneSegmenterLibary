% %load testing workspace
% clear all;
% load('.\seperation_parameters_workspace.mat');
% 
% %using x2
% gf1 = imgaussfilt3(ct3D,x2(1) );
% gf2 = imgaussfilt3(ct3D,x2(2));
% dog_filt = gf1-gf2;
% 
% gf3 = imgaussfilt3(ct3D,x2(3) );
% gf4 = imgaussfilt3(ct3D,x2(4));
% dog_filt2 = gf3-gf4;
% 
% dgcomb = dog_filt+dog_filt2;
% 
% %try nearby_max then comparison
% %using sub-sampling
dg_comb_near_bone = dog_filt;
dg_comb_near_bone(ct3D<1000) = nan;

%nn filling
%KNN FILL NANS
        %knn filling blanks
        currently_labelled = zeros(size(dg_comb_near_bone));
        

        currently_labelled(~isnan(dg_comb_near_bone)) = 1;
        unlabelled_map = zeros(size(dg_comb_near_bone));
        unlabelled_map(isnan(dg_comb_near_bone)) = 1;


        %find nearest point in labelled map for each point in
        %unlabelled_map
        unlabelled_point_indices = find(unlabelled_map);
        [ulp_1,ulp_2,ulp_3] = ind2sub(size(unlabelled_map),unlabelled_point_indices);

        labelled_point_indices = find(currently_labelled);
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
        lookup_n = dg_comb_near_bone(labelled_point_indices);

        dg_comb_near_bone(unlabelled_point_indices) = lookup_n(IDX);



s1=15;
s2=15;
s3=10;
d1_indices = round(1:s1:512-s1);
d2_indices = round(1:s2:512-s2);
d3_indices = round(1:s3:361-s3);
[mg1,mg2,mg3] = meshgrid(d1_indices+(s1/2),d2_indices+(s2/2),d3_indices+(s3/2));
mgI = sub2ind(size(ct3D),mg1,mg2,mg3);
ref_vals_map = nan(size(mgI));
for i_point = 1:length(mgI(:))
    tic
    [mgI_sub1,mgI_sub2,mgI_sub3] = ind2sub(size(mgI),i_point);
    d1min=round(max(1,d1_indices(mgI_sub1)-(0*s1)));
    d1max=round(min(512,d1_indices(mgI_sub1)+(1*s1)));
    d2min=round(max(1,d2_indices(mgI_sub2)-(0*s2)));
    d2max=round(min(512,d2_indices(mgI_sub2)+(1*s2)));
    d3min=round(max(1,d3_indices(mgI_sub3)-(0*s3)));
    d3max=round(min(361,d3_indices(mgI_sub3)+(1*s3)));
    sub_region = dg_comb_near_bone(d1min:d1max,d2min:d2max,d3min:d3max);
    %maxes(i_point) = min(sub_region(:));
    prct(i_point,:) = prctile(sub_region(:),[95 85 75]);
    toc
end

 ref_vals_map(:) = prct(:,1);
% %interp out
[mg1full,mg2full,mg3full]  = meshgrid(1:512,1:512,1:361);
GI = griddedInterpolant({d1_indices+(s1/2),d2_indices+(s2/2),d3_indices+(s3/2)},ref_vals_map, 'linear');
 full_ref_vals = GI({1:512,1:512,1:361});
% %relative_nearby_max = dgcomb-full_maxes;
figure('Name',"prctile_95");
imshow3D(full_ref_vals)

%%% Testing locations
%1.Example rib (left side) slice 0
%2.sternum slice 0, top left
%3.sternum slice 0, bottom right
%4. sternum slice 24, left
%[location_number,info]
%info:
%1: slice
%2: d1_min
%3: d1_max
%4: d2_min
%5: d2_max
testing_locations(1,:) = [1 72 148 106 196];
testing_locations(2,:) = [1 216 119 234 138];
testing_locations(3,:) = [1 258 111 273 128];
testing_locations(4,:) = [24 215 95 232 121];
testing_locations(5,:) = [22 127 109 149 124];
testing_locations(6,:) = [1 389 123 418 164];
for i_tl = 1:size(testing_locations,1)
    cur_tl = testing_locations(i_tl,:);
    %show dgcomb
    figure('Name',"base_edge_map"+num2str(i_tl))
    segment = dog_filt(cur_tl(3):cur_tl(5),cur_tl(2):cur_tl(4),cur_tl(1));
    imshow(segment,[min(segment(:)),max(segment(:))],'InitialMagnification',800,'Border','tight');
    %show prctiles
    for i_pct = 1:3
        ref_vals_map(:) = prct(:,i_pct);
        %interp out
        [mg1full,mg2full,mg3full]  = meshgrid(1:512,1:512,1:361);
        full_ref_vals = interp3(mg1,mg2,mg3,ref_vals_map,mg1full,mg2full,mg3full, 'linear');
        f = figure('Name',"prctile_i_" + num2str(i_pct));
        segment = full_ref_vals(cur_tl(3):cur_tl(5),cur_tl(2):cur_tl(4),cur_tl(1));
        imshow(segment,[min(segment(:)),max(segment(:))],'InitialMagnification',800,'Border','tight');
        waitfor(f)
    end
end


%results (for 70 percentile)
dg_thresh_y = [-750 -325 -300 -240 -500 ...
     -700 -700 -700 -700 ...
     0 0 0 0 0];
prctile_70_x = [250 200 6 -26 150 ...
     400 500 600 700 ...
     -100 -200 -300 -400 -500];
T = table(prctile_70_x',dg_thresh_y');





ref_vals_map = nan(size(mgI));
ref_vals_map(:) = prct(:,10);

[mg1full,mg2full,mg3full]  = meshgrid(1:512,1:512,1:361);
 full_ref_vals = interp3(mg1,mg2,mg3,ref_vals_map,mg1full,mg2full,mg3full, 'linear');
dg_thresh_vals = nan(size(ct3D));
dg_thresh_vals(:) = trainedModel.predictFcn(table(full_ref_vals(:)));

surfaces_map = dgcomb>dg_thresh_vals;
surf_map_dilated = imfilter(mod((double(surfaces_map)+1),2),[1,1,1,1,1;1,1,1,1,1;1,1,1,1,1]);




%% strat2:
%from good dog (x = [0.0104 1.9997])
%threshold ~ 50

surfaces_map = dog_filt>0;
surfaces_map(bone_model_map==0) = 0;
CC = bwconncomp(surfaces_map);
for i_cc = 1:length(CC.PixelIdxList)
    if(length(CC.PixelIdxList{i_cc})<=5)
        surfaces_map(CC.PixelIdxList{i_cc}) = 0;
    end
end

%slice by slice: if not internal then posible dilate location
%2 rounds of r=1 dilation
for i_slice = 1:size(ct3D,3)
    tic
    filled_slice = imfill(surfaces_map(:,:,i_slice),'holes');
    dil = imdilate(surfaces_map(:,:,i_slice),[1 1 1; 1 1 1]);
    dil(filled_slice==1) = 0;
    surfaces_map(:,:,i_slice) = surfaces_map(:,:,i_slice) + dil;
    %round 2
    filled_slice = imfill(surfaces_map(:,:,i_slice),'holes');
    dil = imdilate(surfaces_map(:,:,i_slice),[1 1 1; 1 1 1]);
    dil(filled_slice==1) = 0;
    surfaces_map(:,:,i_slice) = surfaces_map(:,:,i_slice) + dil;
    toc
end

surf_map_dilated = imfilter(surfaces_map,[1,1,1,1,1;1,1,1,1,1;1,1,1,1,1]);

smg3 = surf_map_dilated>=3;
