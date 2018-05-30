%
torso_datapath = Configuration.ExampleDicomDataPath;
head_datapath = Configuration.ExampleHeadDataPath;

[Slices, data3D] = LoadDICOM(torso_datapath, true);



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

hist_bins = [0:1:largest_value];

%[hist1] = histogram(Slices(1).Data, hist_bins);

%imhist_arg = ceil(largest_value * bin_locs(2));

overall_intensity_distribution = zeros(1,largest_value);

individual_histograms = cell(length(Slices),1);

for i = 1:length(Slices)
    %individual_histograms{i} = imhist(Slices(i).Data, imhist_arg);
    cur_hist = histogram(Slices(i).Data, hist_bins);
    individual_histograms{i} = cur_hist.BinCounts;
    overall_intensity_distribution = overall_intensity_distribution + individual_histograms{i};
end


%% Find major peaks in overall distribution
figure
hold on
for i=5:5:25
    plot(individual_histograms{20+i});
end


%super smooth the overall intensity distribution
smoothed_overall_in_dist = smoothdata(overall_intensity_distribution);
smoothed_slice_in_distr = {};
for i = 1:length(Slices)
    smoothed_slice_in_distr{i} = smoothdata(individual_histograms{i}, 'gaussian', 40);
    
    %find first min
    [pks,locs,w,p] = findpeaks(-1*smoothed_slice_in_distr{i}, 1:largest_value);
    
    %find first minima from left and crop the data to remove the large
    %near-zero half-peak
    
    smoothed_slice_in_distr{i}(1:locs(1)) = nan;
end

%find first minima from left and crop the data to remove the large
%near-zero half-peak

%Re-smooth the data to contain only two peaks (from the two groups of
%peaks)

%locate the two peaks for each slice
plot_figuresA = false;
if (plot_figuresA)
figure
plot(smoothed_overall_in_dist);
title('overall');
xlabel('CT data intensity value');
ylabel('Quantity of CT datapoints (in this slice, at specified intensity value)');

figure
plot(smoothed_slice_in_distr{200});
title('Slice 200');
xlabel('CT data intensity value');
ylabel('Quantity of CT datapoints (in this slice, at specified intensity value)');

figure
plot(smoothed_slice_in_distr{190});
title('Slice 190');
xlabel('CT data intensity value');
ylabel('Quantity of CT datapoints (in this slice, at specified intensity value)');

figure
plot(smoothed_slice_in_distr{210});
title('Slice 210');
xlabel('CT data intensity value');
ylabel('Quantity of CT datapoints (in this slice, at specified intensity value)');

figure
plot(smoothed_slice_in_distr{180});
title('Slice 180');
xlabel('CT data intensity value');
ylabel('Quantity of CT datapoints (in this slice, at specified intensity value)');

figure
plot(smoothed_slice_in_distr{170});
title('Slice 170');
xlabel('CT data intensity value');
ylabel('Quantity of CT datapoints (in this slice, at specified intensity value)');

figure
plot(smoothed_slice_in_distr{195});
title('Slice 195');
xlabel('CT data intensity value');
ylabel('Quantity of CT datapoints (in this slice, at specified intensity value)');
end
%% Peak Identification

%p1-5 as described in G_0__RoughOutline.odt

[   peak_locations,...
    coordinate_transform_parameters,...
    overall_distr_peak_locations,...
    overall_distr_transform_parameters] = IdentifyPeaks(smoothed_overall_in_dist, smoothed_slice_in_distr);

