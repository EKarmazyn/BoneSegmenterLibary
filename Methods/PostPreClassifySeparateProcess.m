function [connected_volumes_array] = PostPreClassifySeparateProcess(final_labelled_map,taggedIdx_array,ct3D,dataset,X_ct,SCALE )
%POSTPRECLASSIFYSEPARATEPROCESS extra processing before classification

 zero_map = uint8(zeros(size(ct3D)));
 zero_map16 =  uint16(zeros(size(ct3D)));

initial_num_labels = max(final_labelled_map(:));

%new_label_counter = 1;


%generate old_pixelidxlist
%old_pixelidxlist= cell(initial_num_labels,1);


all_indices = cell2mat(taggedIdx_array');

indice_labels = final_labelled_map(all_indices);


% for i_label = 1:initial_num_labels
%     old_pixelidxlist{i_label} = all_indices(indice_labels==i_label);
% end

old_pixelidxlist = taggedIdx_array;

full_pixelIdx_list = cell(initial_num_labels,1);
bound_props = cell(initial_num_labels,1);
parfor i_label = 1:initial_num_labels
    tic
    cur_map = zero_map;
    cur_map(old_pixelidxlist{i_label}) = 1;
    [bound_map,ind1,ind2,ind3, full_size] = BoundMap(cur_map,0);
    %ind are the offsets so need to ind2sub to convert
    CC_cur = bwconncomp(bound_map);
    
    if(length(CC_cur.PixelIdxList)>1)
        %adjust labels
        %[s1, s2, s3] = ind2sub(CC_cur.PixelIdxList,size(bound_map));
        %full_pixelIdx_list{i_label} = sub2ind(size(cur_map),s1+ind1,s2+ind2,s3+ind3);
        bound_props{i_label}.Inds = [(ind1-1),(ind2-1),(ind3-1)];
        bound_props{i_label}.Size = size(bound_map);
        full_pixelIdx_list{i_label} = CC_cur.PixelIdxList;
        
    elseif(length(CC_cur.PixelIdxList)==1)
        [s1, s2, s3] = ind2sub(size(bound_map),CC_cur.PixelIdxList{1});
        full_pixelIdx_list{i_label} = sub2ind(size(cur_map),s1+(ind1-1),s2+(ind2-1),s3+(ind3-1));
       
        %full_pixelIdx_list{i_label} = CC_cur.PixelIdxList{1};
    else
        full_pixelIdx_list{i_label} = [];
    end
    toc
end

full_long_pixelIdx_list = cell(initial_num_labels*20,1);


new_label_counter = 1;
%update final_labelled_map
for i_label = 1:initial_num_labels
    if(iscell(full_pixelIdx_list{i_label}))
        %need to update
        [s1, s2, s3] = ind2sub(bound_props{i_label}.Size,full_pixelIdx_list{i_label}{1});
        inds = bound_props{i_label}.Inds;
        full_long_pixelIdx_list{i_label} = sub2ind(size(zero_map),s1+inds(1),s2+inds(2),s3+inds(3));
       
        
        
        %full_long_pixelIdx_list{i_label} = full_pixelIdx_list{i_label}{1};
        for i_new_label = 2:length(full_pixelIdx_list{i_label})
            
            [s1, s2, s3] = ind2sub(bound_props{i_label}.Size,full_pixelIdx_list{i_label}{i_new_label});
            inds = bound_props{i_label}.Inds;
            full_long_pixelIdx_list{initial_num_labels+new_label_counter} = sub2ind(size(zero_map),s1+inds(1),s2+inds(2),s3+inds(3));
       
            
            
            final_labelled_map(full_long_pixelIdx_list{initial_num_labels+new_label_counter}) = initial_num_labels+new_label_counter;
            %full_long_pixelIdx_list{initial_num_labels+new_label_counter} = full_pixelIdx_list{i_label}{i_new_label};
        
            new_label_counter = new_label_counter+1;
        end
    else
        %dont need to update
        full_long_pixelIdx_list{i_label} = full_pixelIdx_list{i_label};
        
    end
end



new_num_vols = max(final_labelled_map(:));
%double check:
d_check = length(unique(final_labelled_map(:))) -1;

full_long_pixelIdx_list = full_long_pixelIdx_list(~cellfun(@isempty, full_long_pixelIdx_list));





%split volumes into large, medium and small
full_long_pixelIdx_list = full_long_pixelIdx_list(~cellfun(@isempty, full_long_pixelIdx_list));

vol_lengths = zeros(new_num_vols,1);

%of the small vols,
%   find those without other 'bone' nearby and remove
for i_label = 1:new_num_vols
    vol_lengths(i_label) = length(full_long_pixelIdx_list{i_label});
end

short_vol_min = 237.8906*0.3965; %100 for 1

short_label_min = round(short_vol_min/SCALE.vol);
short_labels = find(vol_lengths<short_label_min);

 rp  = regionprops3(final_labelled_map, {'Centroid','EquivDiameter'});
 CoM = rp(:,'Centroid');
 EqDia = rp(:,'EquivDiameter');
%setup co-ords transform
%gi = griddedInterpolant(1:size(ct3D,1), 1:size(ct3D,2), 1:size(ct3D,3), ndgrid(X_ct{1},X_ct{2},X_ct{3}));
gi = {};
for i = 1:3
    gi{i} = griddedInterpolant(1:size(ct3D,i), X_ct{i});

end

all_indices = cell2mat(full_long_pixelIdx_list);


small_indices_arr = {};

for i_short = 1:length(short_labels)
    cur_label = short_labels(i_short);
        small_indices_arr{i_short} = full_long_pixelIdx_list{cur_label};

end

%remove cur label indices from all indices list
all_small_indices = cell2mat(small_indices_arr');
other_indices = setdiff(all_indices, all_small_indices);

%use NN tagging to find distance to closest part of the model

%unlabelled_point_indices = setdiff(unlabelled_point_indices,all_seed_indexes);

[ul_arr] = CoM{short_labels, 'Centroid'};

%labelled_point_indices = find(labelled_map);
%labelled_point_indices = all_seed_indexes;
[lp_1,lp_2,lp_3] = ind2sub(size(final_labelled_map),other_indices);
%convert to real positions using x_ct
ulp_1 = gi{1}(ul_arr(:,2));%swapped on purpose!
ulp_2 = gi{2}(ul_arr(:,1));%
ulp_3 = gi{3}(ul_arr(:,3));

lp_1 = X_ct{1}(lp_1);
lp_2 = X_ct{2}(lp_2);
lp_3 = X_ct{3}(lp_3);



[~, D] = knnsearch([lp_1',lp_2',lp_3'],[ulp_1,ulp_2,ulp_3]);


%f = figure
binary_labelled_map = zero_map;
binary_labelled_map(all_indices) = 1;

short_rs =  EqDia{short_labels, 'EquivDiameter'};
short_rs = short_rs/2;

D_min_r = D-short_rs;

%dist_range = 1:80;
min_allowed_distance = 3;
dist_inds = find(D_min_r>=min_allowed_distance); %to_purge
%dist_inds = find(D_min_r<3); % survigin
culled_labels = short_labels(dist_inds);
% target_vol = sel_labels;
% use_indices = true;
% target_vol_inds = {};
% for i = 1:length(sel_labels)
%     target_vol_inds{i} = full_long_pixelIdx_list{sel_labels(i)};
% end


%view_individual_volume_script;
valid_volumes = setdiff(1:length(full_long_pixelIdx_list),culled_labels);


%split volumes into large, medium and small for easier classification!
valid_lengths = vol_lengths(valid_volumes);

v_small_boundary = round(10*(2.3789*0.3965/SCALE.vol));
v_small_vols = valid_volumes(find(valid_lengths<=v_small_boundary));

small_boundary = round(50*(2.3789*0.3965/SCALE.vol));
small_vols = valid_volumes(find(valid_lengths<=small_boundary));
only_small_vols = setdiff(small_vols, v_small_vols);




medium_boundary = round(200*(2.3789*0.3965/SCALE.vol));
medium_vols = valid_volumes(find(valid_lengths<=medium_boundary));
only_medium_vols = setdiff(medium_vols, small_vols);

large_boundary = round(500*(2.3789*0.3965/SCALE.vol));
large_vols = valid_volumes(find(valid_lengths<=large_boundary));
only_large_vols = setdiff(large_vols, medium_vols);


v_large_vols = valid_volumes(find(valid_lengths>large_boundary));
only_v_large_vols = setdiff(v_large_vols, large_vols);

%make a label map by vol size for good picture
size_label_map = zero_map;


cur_list = v_small_vols;
cur_lab = 1;
for i_v = 1:length(cur_list)
    size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
end

cur_list = only_small_vols;
cur_lab = 2;
for i_v = 1:length(cur_list)
    size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
end

cur_list = only_medium_vols;
cur_lab = 3;
for i_v = 1:length(cur_list)
    size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
end

cur_list = only_large_vols;
cur_lab = 4;
for i_v = 1:length(cur_list)
    size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
end

cur_list = only_v_large_vols;
cur_lab = 5;
for i_v = 1:length(cur_list)
    size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
end

cur_list = culled_labels;
cur_lab = 6;
for i_v = 1:length(cur_list)
    size_label_map(full_long_pixelIdx_list{cur_list(i_v)}) = cur_lab;
end

% f% = figure
% [hpat] = ViewSeperateVolumes(size_label_map, X_ct);
% 
% for i_h = 1:6
%     hpat{i_h}.FaceAlpha = 0.1;
%     hpat{i_h}.EdgeAlpha = 0.1;
% end

%save("C:\Users\mazna\Documents\nl\U\P\Report\FinalReport\FigureData\HardSeperate\final_" + num2str(dataset),'f');

%close(f);
clear connected_volumes_array;
connected_volumes_array(1,length(valid_volumes)) = ConnectedVolume;

for i_valid_vol = 1:length(valid_volumes)
    connected_volumes_array(1,i_valid_vol) = ConnectedVolume(full_long_pixelIdx_list{valid_volumes(i_valid_vol)});
    
end

%add some already calculated features/properties


cur_list = v_small_vols;
cur_lab = 1;
for i_v = 1:length(cur_list)
    connected_volumes_array(find(valid_volumes==cur_list(i_v))).VolumeSize = cur_lab;
end

cur_list = only_small_vols;
cur_lab = 2;
for i_v = 1:length(cur_list)
    connected_volumes_array(find(valid_volumes==cur_list(i_v))).VolumeSize = cur_lab;
end

cur_list = only_medium_vols;
cur_lab = 3;
for i_v = 1:length(cur_list)
    connected_volumes_array(find(valid_volumes==cur_list(i_v))).VolumeSize = cur_lab;
end

cur_list = only_large_vols;
cur_lab = 4;
for i_v = 1:length(cur_list)
    connected_volumes_array(find(valid_volumes==cur_list(i_v))).VolumeSize = cur_lab;
end

cur_list = only_v_large_vols;
cur_lab = 5;
for i_v = 1:length(cur_list)
    connected_volumes_array(find(valid_volumes==cur_list(i_v))).VolumeSize = cur_lab;
end




end

