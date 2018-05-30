
%split labels using seeds
new_label_idx_array = {};
%cell_ind_counter = 1;
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
        new_label_idx_array{i_cc_sep} = taggedIdx_array;
        %cell_ind_counter = cell_ind_counter + length(taggedIdx_array);
        %toc(TT);
    else
        new_label_idx_array{i_cc_sep} = CC_seperated.PixelIdxList{i_cc_sep};
        %cell_ind_counter = cell_ind_counter+1;
    end
    %toc(T);
end