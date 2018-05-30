function [prob_bone_seed_map, prob_minimum_bone_seed_intensity, def_bone_seed_map, def_bone_seed_min_intensity] = FindBoneSeeds(ct3D, valid_areas_map)
%FINDBONESEEDS Summary of this function goes here
%   should work on 2d or 3d

prob_bone_seed_map = zeros(size(ct3D));
def_bone_seed_map = zeros(size(ct3D));
%bone_seed_value_map = nan(size(ct3D));

%% Overall Intensity Distribution
largest_value = max(max(max(ct3D)));
overall_distribution_obj = histogram(ct3D, [0:1:largest_value]);

overall_distribution = overall_distribution_obj.BinCounts;
close
peak_locations = PrimitiveFindPeaksInDistr(smooth(overall_distribution), false);

%after the last found peak look at the distribution and then calculate the
%bone cut-off
peak_locations(isnan(peak_locations)) = [];
cutoff = peak_locations(end);
dist = overall_distribution(cutoff:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACTUALLY WORK IT OUT?
%phat = mle(dist, 'distribution','hn');
flipdist = flip(dist);
scale = max(dist);
    
cdf = strcat(num2str(scale/2),'*(1+erf((x-a)/(b*sqrt(2))))');
mean_guess = length(dist);
startPoints = [mean_guess 1];
[f1, gof, output] = fit((1:length(flipdist))',flipdist', cdf, 'Start' , startPoints);


%     figure
%     plot(f1)%,(1:dist(dist))',dist')
%     hold on;
%     plot(flipdist);
%     close


prob_var_multiplier = 20;
def_var_multiplier = 60;

bone_seed_threshold = 1600; %OVERRIDE%%%round(cutoff + (f1.b*prob_var_multiplier));
prob_bone_seed_map(ct3D > bone_seed_threshold) = 1;

prob_bone_seed_map(valid_areas_map == 0) = 0;

prob_minimum_bone_seed_intensity = bone_seed_threshold;

bone_seed_threshold = 2300;%OVERRIDE%%%round(cutoff + (f1.b*def_var_multiplier));
def_bone_seed_map(ct3D > bone_seed_threshold) = 1;
def_bone_seed_map(valid_areas_map == 0) = 0;
def_bone_seed_min_intensity = bone_seed_threshold;

end

