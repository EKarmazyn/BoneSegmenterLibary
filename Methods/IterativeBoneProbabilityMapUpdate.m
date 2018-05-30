function [intensity_map, bone_probability_map, fixed_map] = IterativeBoneProbabilityMapUpdate(intensity_map,bone_probability_map, fixed_map)
%ITERATIVEBONEPROBABILITYMAPUPDATE An iteration to update the
%bone_probability_map

%find all values adjacent to previously calculated values
to_calculate_binary_map = false(size(intensity_map));
to_calculate_binary_map(bone_probability_map == 1) = true;
to_calculate_binary_map = imdilate(to_calculate_binary_map, strel('square',3));
to_calculate_binary_map(bone_probability_map == 1) = false;
to_calculate_binary_map(fixed_map == true) = false;





end

