function [final_labelled_map, taggedIdx_array] = PreClassifySeperate2(ct3D,binary_map_filled, CT_dimension_spacing,X_ct,SCALE, use_coating, coating_radius_addition)
%PRECLASSIFYSEPERATE2 Summary of this function goes here
%coating_radius_addition : on top of half of the extra seed avoidance
%radius
z_thresh = 25
xy_smoothed_thresh = 30
[voided, voided_eroded, voided_eroded2] = HardSeperateVolumes(ct3D,binary_map_filled, CT_dimension_spacing, xy_smoothed_thresh,z_thresh, SCALE);


comb = voided;
comb(voided_eroded==1) = 2;
comb(voided_eroded2==1) = 3;
%     
% figure
% imshow3D(comb);


%MULTILEVEL SEPERATION!

%first tag upper level
CC_primary = bwconncomp(voided_eroded2);
labelled_prime_seed_map = labelmatrix(CC_primary);

%secondary level
CC2 = bwconncomp(voided_eroded);
labelled_sec_seed_map = labelmatrix(CC2);

%tertiary level
CC3 = bwconncomp(voided);
labelled_tert_seed_map = labelmatrix(CC3);


%get large primary seeds
small_map = zeros(size(ct3D));
min_primary_length = 70;
for i_CC = 1:length(CC_primary.PixelIdxList)
    if(length(CC_primary.PixelIdxList{i_CC})<min_primary_length)
        small_map(CC_primary.PixelIdxList{i_CC}) = i_CC;
        labelled_prime_seed_map(CC_primary.PixelIdxList{i_CC}) = 0;
    end
    
end
%re-label primary seeds
binary_primary_seeds = zeros(size(ct3D));
binary_primary_seeds(find(labelled_prime_seed_map)) = 1;
CC_primary = bwconncomp(binary_primary_seeds);
labelled_prime_seed_map = labelmatrix(CC_primary);




%for each secondary seed:
valid_secondary_seeds = [];
min_secondary_length = 0;
for i_CC = 1:length(CC2.PixelIdxList)
    %if it connects more than one primary seed then its bad
    %if it contains a primary seed then its a potential additional error
    %and the primary seed should be used instead
    %min length too
    if(length(unique(labelled_prime_seed_map(CC2.PixelIdxList{i_CC})))==1 && ...
        length(CC2.PixelIdxList{i_CC})>     min_secondary_length)
        valid_secondary_seeds = [valid_secondary_seeds; i_CC];
    end
    
end

max_primary_seed = max(labelled_prime_seed_map(:));
valid_secondary_seed_map = uint16(zeros(size(ct3D)));
for i = 1:length(valid_secondary_seeds)
    valid_secondary_seed_map(CC2.PixelIdxList{valid_secondary_seeds(i)}) = max_primary_seed+i;
end

prim_sec_seed_map = uint16(labelled_prime_seed_map) + valid_secondary_seed_map;

%figure
%ViewSeperateVolumes(valid_secondary_seed_map,X_ct);


%repeat with teriary
%for each tertiary seed:
valid_tertiary_seeds = [];
min_tertiary_length = 10;
for i_CC = 1:length(CC3.PixelIdxList)
    %if it connects more than one primary/sec seed then its bad
    %if it contains a single primary/sec seed then its a potential additional error
    %and the primary/sec seed should be used instead
    %min length too
    if(length(unique(prim_sec_seed_map(CC3.PixelIdxList{i_CC})))==1 && ...
        length(CC3.PixelIdxList{i_CC})>     min_tertiary_length)
        valid_tertiary_seeds = [valid_tertiary_seeds; i_CC];
    end
    
end

max_prim_sec_seed = max(prim_sec_seed_map(:));
valid_tertiary_seed_map = uint16(zeros(size(ct3D)));
for i = 1:length(valid_tertiary_seeds)
    valid_tertiary_seed_map(CC3.PixelIdxList{valid_tertiary_seeds(i)}) = max_prim_sec_seed+i;
end

full_seed_map = prim_sec_seed_map + valid_tertiary_seed_map;

% figure
% ViewSeperateVolumes(prim_sec_seed_map,X_ct);

%use the super hard seperate cvs to connect generated seeds, t


%super_voided = HardSeperateVolumesOld(ct3D,binary_map_filled, CT_dimension_spacing, 15,0);

%%% SS2 seperate into distinct volumes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%generate edge map
dog_gf1 = 0.0104;
dog_gf2 = 1.9997;
dog_filt = imgaussfilt3(double(ct3D),dog_gf1)-imgaussfilt3(double(ct3D),dog_gf2);
%tic
smoothed = imgaussfilt3(dog_filt,4,'FilterSize',11);
%toc

dog_minus_smoothed = dog_filt - smoothed;

basic_surfaces = dog_minus_smoothed>50;


%
updated_map = binary_map_filled;
map_coating = imdilate(updated_map,[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1]);
map_coating(updated_map==1) = 0;



%map_plus_coating = updated_map;
%map_plus_coating = uint8(map_plus_coating);
%map_plus_coating(map_coating==1) = 1;


%seperate map plus coating
seperated_map = updated_map;
seperated_map(basic_surfaces==1) = 0;

eroded_seperated_map = imerode(seperated_map, strel('sphere',2));


%figure
%ViewSeperateVolumes(labelled_seperated, X_ct);


%method to add some extra seed labels in missed regions, (particularly
%ilium and sections of rib)

%heavily erode the filled_bone_binary_map which has clear ribs/illium 
heavily_eroded_binary_map = imerode(binary_map_filled, strel('sphere',2));
% figure
% ViewSeperateVolumes(heavily_eroded_binary_map, X_ct);

%remove points near a dilated current_label_map
%need coating to be further away i.e. half of this dilation needs to be
%less than the coating layer 1 radius
labelled_seperated_binary = eroded_seperated_map;
labelled_seperated_binary(eroded_seperated_map~=0) = 1;
r1 = 6;
dilated_labelled = imdilate(labelled_seperated_binary, strel('sphere',r1));
extra_seeds = heavily_eroded_binary_map;
extra_seeds(dilated_labelled==1) = 0;

%figure('Name','extra seeds');
%ViewSeperateVolumes(extra_seeds,   X_ct);

%set as seed_labels
eroded_seperated_map(extra_seeds == 1) = 1;


CC_seperated = bwconncomp(eroded_seperated_map);
labelled_seperated=labelmatrix(CC_seperated);


%split labels using seeds
new_label_idx_array = {};
for i_cc_sep = 1:length(CC_seperated.PixelIdxList)
    %T = tic;
    %if it contains multiple seeds it needs to be split. 
    inc_seeds = unique(full_seed_map(CC_seperated.PixelIdxList{i_cc_sep}));
    
    if(length(inc_seeds)>2)% 2 is zero + seed label
        % seed map section
        inc_seed_array = {};
        for i_inc_seed = 2:length(inc_seeds) % 1 is 0
            inc_seed_array{i_inc_seed-1} = CC_seperated.PixelIdxList{inc_seeds(i_inc_seed)};
        end
        %TT = tic;
        taggedIdx_array = NNTaggingFromIndices(inc_seed_array, CC_seperated.PixelIdxList{i_cc_sep}, X_ct, size(ct3D));
        new_label_idx_array{i_cc_sep} =  taggedIdx_array;
        %toc(TT);
    else
        new_label_idx_array{i_cc_sep} = CC_seperated.PixelIdxList{i_cc_sep};
    end
    %toc(T);
end



if(use_coating)
%%%%%%%%%%%%%%-----------------
%errors can occur when there is no label in a non bone region which is
%included in the binary bone map, therefore coating needs to be used on the
%seperated map !, if used on the eroded seperated then coating is added to
%bone regions where there is thin bone,
r2 = r1/2+coating_radius_addition;
under_coat = seperated_map + extra_seeds;
coating_layer_1 = imdilate(under_coat, strel('sphere',r2));
coating_layer_2 = imdilate(coating_layer_1, strel('sphere',3)) - coating_layer_1;

surfaces_avoiding = imdilate(basic_surfaces, strel('sphere',2));

coating_final = coating_layer_2;
coating_final(surfaces_avoiding==1) = 0;
%figure('Name','coating_final') 
%ViewSeperateVolumes(coating_final, X_ct);
end

%convert new_label_idx_array into long cell array
%COULD POTENTIALLY FILTER OUT SMALL HERE
long_label_idx_array = {};
cell_ind = 1;
for i = 1:length(new_label_idx_array)
    if(iscell(new_label_idx_array{i}))
        for j = 1:length(new_label_idx_array{i})
            long_label_idx_array{cell_ind} = new_label_idx_array{i}{j};
            cell_ind = cell_ind+1;
        end
    else
        long_label_idx_array{cell_ind} = new_label_idx_array{i};
        cell_ind = cell_ind+1;
    end
end

if(use_coating)
%label coating
coating_label = length(long_label_idx_array) + 1;
long_label_idx_array{coating_label} = find(coating_final);
end


%KNN fill eroded points using labels
taggedIdx_array = NNTaggingFromIndices(long_label_idx_array, find(eroded_seperated_map-seperated_map), X_ct, size(ct3D));


cur_labelled_map = uint16(zeros(size(ct3D)));
for i_label = 1:length(taggedIdx_array)
    cur_labelled_map(taggedIdx_array{i_label}) = i_label;
end

%cur_labelled_map(cur_labelled_map==coating_label) = 0;
% figure('Name','cur_labelled_map_minus_coating') 
% ViewSeperateVolumes(cur_labelled_map, X_ct);


%KNN fill hard bone points  using labels
min_hard_bone = 1500;
hard_bone_idx = find(ct3D>min_hard_bone);
all_labelled_ind = cell2mat(taggedIdx_array');
hard_bone_idx = setdiff(hard_bone_idx, all_labelled_ind);
taggedIdx_array = NNTaggingFromIndices(taggedIdx_array, hard_bone_idx, X_ct, size(ct3D));


cur_labelled_map = uint16(zeros(size(ct3D)));
for i_label = 1:length(taggedIdx_array)
    cur_labelled_map(taggedIdx_array{i_label}) = i_label;
end

%view_map = cur_labelled_map;
%view_map(cur_labelled_map==coating_label) = 0;

% figure('Name','after_hard_bone_tagging_labelled_map_minus_coating') 
% ViewSeperateVolumes(view_map, X_ct);
% 

%KNN fill rest of the seperated dataset points  using labels
labelled_seperated(cur_labelled_map~=0) = 0;
lab_sep_idc = find(labelled_seperated);
taggedIdx_array = NNTaggingFromIndices(taggedIdx_array, lab_sep_idc, X_ct, size(ct3D));



cur_labelled_map = uint16(zeros(size(ct3D)));
for i_label = 1:length(taggedIdx_array)
    cur_labelled_map(taggedIdx_array{i_label}) = i_label;
end

%view_map = cur_labelled_map;
%view_map(cur_labelled_map==coating_label) = 0;
%
% figure('Name','after_seperated_filled_minus_coating') 
% ViewSeperateVolumes(view_map, X_ct);

%some sections have been incorrectly markes as coating, so we need
%additional beginning labels! DONE!

%need to carefully re-map definate bone points away from coating
%no longer needed  ( i think!)


%potential thinning step here? all the bone basically has a certain density
%unhelpful at this stage as it encourages coating errors


%check coating is doing well:
% view_slice = squeeze(cur_labelled_map(:,200,:));
% view_slice(view_slice~= 0) = 1;
% figure
% imshow(view_slice, [0 1]);

%KNN fill rest of the binary_bone_map using labels
binary_bone_map_indixes = find(binary_map_filled);

all_labelled_ind = cell2mat(taggedIdx_array');
binary_bone_map_indixes = setdiff(binary_bone_map_indixes, all_labelled_ind);
taggedIdx_array = NNTaggingFromIndices(taggedIdx_array, binary_bone_map_indixes, X_ct, size(ct3D));

cur_labelled_map = uint16(zeros(size(ct3D)));
for i_label = 1:length(taggedIdx_array)
    cur_labelled_map(taggedIdx_array{i_label}) = i_label;
end



final_labelled_map = cur_labelled_map;
if(use_coating)
    final_labelled_map(cur_labelled_map==coating_label) = 0;
    taggedIdx_array{coating_label} = [];
end

return
% figure('Name','final_seperated') 
% ViewSeperateVolumes(final_labelled_map, X_ct);



end

