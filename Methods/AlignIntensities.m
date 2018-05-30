function [CorrectedSlices, AlignmentCoefficients] = AlignIntensities(Slices, data3D)
%ALIGNINTENSITIES Aligns the intensity distributions of the input slices
%(of CT data assumed)
%   Slices  : Vector of Slice structures which have Info, Data and
%   Position. 
%   Alignment Coefficients : A matrix (slice, coefficent) containing:
%       1: Position of peak region 1 (in normalised co-ordinates)
%       2: Position of peak region 2 (in normalised co-ordinates)
%
%
%   Normalised co-ordinate system : 
%                                   0 -> 0, 
%                                   1 -> first peak from whole dataset, 
%                                   2 -> second peak from whole dataset



%Clone slices for corrected
CorrectedSlices = [];
for i = 1:length(Slices)
    CorrectedSlices(1,i).Info = Slices(1,i).Info;
    CorrectedSlices(1,i).Data = Slices(1,i).Data;
    CorrectedSlices(1,i).Position = Slices(1,i).Position;
end


%Alignment coefficents
AlignmentCoefficients = nan(length(Slices),2);

largest_value = max(max(max(data3D)));

%% Partition overall intensity distribution (histogram) 

[hist1, bin_locs] = imhist(Slices(1).Data, largest_value);

imhist_arg = ceil(largest_value * bin_locs(2));

overall_intensity_distribution = zeros(largest_value,1);

individual_histograms = cell(length(Slices),1);

for i = 1:length(Slices)
    individual_histograms{i} = imhist(Slices(i).Data, imhist_arg);
    individual_histograms{i} = individual_histograms{i}(1:largest_value);
    overall_intensity_distribution = overall_intensity_distribution + individual_histograms{i};
end


%% Find two major peaks in overall distribution

%super smooth the overall intensity distribution
smoothed_overall_in_dist = smoothdata(overall_intensity_distribution);
smoothed_slice_in_distr = {};
for i = 1:length(Slices)
    smoothed_slice_in_distr{i} = smoothdata(individual_histograms{i});
end

%find first minima from left and crop the data to remove the large
%near-zero half-peak

%Re-smooth the data to contain only two peaks (from the two groups of
%peaks)

%locate the two peaks for each slice

figure
plot(smoothed_overall_in_dist);
title('overall');
figure
plot(smoothed_slice_in_distr{1});
title('1');
figure
plot(smoothed_slice_in_distr{10});
title('10');

figure
plot(smoothed_slice_in_distr{100});
title('100');

figure
plot(smoothed_slice_in_distr{400});
title('400');



[overall_peaks, locs] = findpeaks(overall_intensity_distribution, 1:largest_value);

end

