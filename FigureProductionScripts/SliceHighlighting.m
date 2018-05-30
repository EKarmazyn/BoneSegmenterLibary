%[Slices, data3D] = LoadDICOM(Configuration.ExampleDicomDataPath, true);


% 
% %Clone slices for corrected
% CorrectedSlices = [];
% for i = 1:length(Slices)
%     CorrectedSlices(1,i).Info = Slices(1,i).Info;
%     CorrectedSlices(1,i).Data = Slices(1,i).Data;
%     CorrectedSlices(1,i).Position = Slices(1,i).Position;
% end


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

%% Show range in highlighted colour

slice_to_show = 10;
%range = [28 71]; %peak 1 (for slice 1)
%range = [71 200]; %peak 2
%range = [863 974]; %peak 3 
%range = [974 1142]; %peak 4 
range = [1200 largest_value]; %peak 5 
figure

imshow(Slices(slice_to_show).Data, [0 largest_value]);
hold on;
highlighted_area = zeros(size(Slices(slice_to_show).Data));

highlighted_area(Slices(slice_to_show).Data >= range(1) & Slices(slice_to_show).Data < range(2)) = largest_value;

g = imshow(highlighted_area, [0 largest_value]);

transparency_map = ones(size(highlighted_area));
transparency_map(0 == (highlighted_area)) = 0;
hold off;
set(g, 'AlphaData', transparency_map);
title('Peak 5 areas');
