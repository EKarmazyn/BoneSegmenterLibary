function [marked_cv_array, binary_map_filled] = VolumeClassifier(connected_volumes_array,X_ct, sizeImg)
%VOLUMECLASSIFIER Summary of this function goes here



fully_labelled_map = uint16(zeros(sizeImg));
for i = 1:length(connected_volumes_array)
    
     
    
    fully_labelled_map(connected_volumes_array(i).IndVoxelList) = i;
    
        
end



ss_lab_map = uint16(zeros(sizeImg));
cur_size = 5; 
cur_list_counter = 1;
cva_lengths = [];
super_short_list = [];
for i = 1:length(connected_volumes_array)
    cva_lengths = [cva_lengths; connected_volumes_array(i).NumVoxels];
    
    
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




%local solidity features? i.e. how much as % is lost on small erode?
binary_map = full_wo_ss_label;
binary_map(binary_map~=0 ) = 1;

ConnectedVolume.GenerateErodeFraction(connected_volumes_array, binary_map);


%red-green [0 0 1] -> [1 0 0]
for i_va = 1:length(connected_volumes_array)
    frac = connected_volumes_array(i_va).Features.ErodeFraction;
    frac_a(i_va) = frac;
end


erode_valid_ids = find(frac_a<0.7);
final_model = uint16(zeros(sizeImg));
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

SaveBone3D(sizeImg,connected_volumes_array, X_ct);

marked_cv_array = connected_volumes_array;


end

