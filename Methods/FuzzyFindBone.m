function [fuzzy_bone_map, initial_threshold, probable_min_bone_intensity] = FuzzyFindBone(ct3D, valid_areas_map)
%FINDBONESEEDS Summary of this function goes here
%   should work on 2d or 3d

if(nargin<3)
    dilate = false;
    local_adaptive = false;
end

fuzzy_bone_map = zeros(size(ct3D));
fuzzy_bone_map = uint8(fuzzy_bone_map);
%bone_seed_value_map = nan(size(ct3D));

%% Overall Intensity Distribution
largest_value = max(max(max(ct3D)));
overall_distribution_obj = histogram(ct3D, [0:1:largest_value]);
overall_distribution = overall_distribution_obj.BinCounts;
smoothed = smooth(overall_distribution);
[peak_locations, soft_tissue_std] = PrimitiveFindPeaksInDistr(smoothed, false);

%after the last found peak look at the distribution and then calculate the
%bone cut-off
% if(any(isnan(peak_locations)))
%     shift_right_adjust=true;
% 
% else
%     shift_right_adjust=false;
% end
% peak_locations(isnan(peak_locations)) = [];
% cutoff = peak_locations(end);
% dist = overall_distribution(cutoff:end);


initial_threshold = peak_locations(4) + soft_tissue_std*8;
fuzzy_bone_map(ct3D > initial_threshold) = 1;

%probable min bone intensity can be hidden by the soft tissue peak quite a
%lot...
%its often safer to just use a fixed value
probable_min_bone_intensity = 850;

% bone_fuzzy_threshold = round(cutoff + 100);
% if(shift_right_adjust)
%     bone_fuzzy_threshold = bone_fuzzy_threshold+100;
% end
% 
% fuzzy_bone_map(ct3D > bone_fuzzy_threshold) = 1;
% 
% fuzzy_bone_map(valid_areas_map == 0) = 0;
% 
% minimum_bone_intensity = (cutoff + 100);







end

