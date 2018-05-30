function [peak_locations, coordinate_transform_parameters, overall_distr_peak_locations, overall_distr_transform_parameters] = IdentifyPeaks(smoothed_overall_in_dist, smoothed_slice_in_distr)
%IDENTIFYPEAKS Summary of this function goes here
%   smoothed_overall_in_dist :: Smoothed version of the overall intensity
%   distribution (pre-normalisation)
%   smoothed_slice_in_distr  :: Cell array of the smoothed intensity
%   distributions of each slice
%   peak_locations           :: Cell array, one entry for each slice.
%   peak_locations{slice}(peak_num). Gives the location of a peak in
%   co-ordinate system where 1 = peak1/2 CoM, 2 = peak 4.

%% First identify peaks and their locations in the overall in_dist
% This has been tested for: 
% This has not been tested for: HEAD, UPPER TORSO, LOWER TORSO, PELIC, LEGS

display = true;
overall_distr_peak_locations = PrimitiveFindPeaksInDistr(smoothed_overall_in_dist, display);
% calculate transform parameters;

%loc of 1.0: 
overall_distr_transform_parameters(1) = overall_distr_peak_locations(2);
%loc of 2.0
overall_distr_transform_parameters(2) = overall_distr_peak_locations(4);
%range between 1.0 and 2.0
unit_transform = overall_distr_peak_locations(4) - overall_distr_peak_locations(2);

%transform to new co-ordinate system
origin = overall_distr_transform_parameters(2) - 2*unit_transform;
for p = 1:5
    overall_distr_peak_locations(p) = (overall_distr_peak_locations(p)-origin)/unit_transform;
end




%% Indentify peaks in each slice
display = false;

for slice_no = 1:length(smoothed_slice_in_distr)
    
    peak_locations{slice_no} = PrimitiveFindPeaksInDistr(smoothed_slice_in_distr{slice_no}, display, strcat('slice',num2str(slice_no)));
end

%im some cases peak 3 is identified as 4, and therefore p4 as p5. To
%correct for this we can use the fact that we know the slices are adjacent
%to each other, and therefore the intensity distributions of locally close
%slices should be similar.

prev_peak_locations = peak_locations{1};

show_fixes = true;

for slice_no = 2:length(smoothed_slice_in_distr)
    
    for p_num = 1:5
        %if peak appears to have drifted by too much, shift peak
        %identifications to fix
        peak_velocity = (peak_locations{slice_no}(p_num) - prev_peak_locations(p_num))/1; %can divide by inter-slice distance here
        if(peak_velocity > Configuration.PeakAlignment_MaxPeakVelocity)
            %find the peak it should have been labelled as
            %%%%%%%%%%%%%%%%%%%%%%%% HEREHEREHEREHEREHERE
            %%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        end
        %if fixed, show old and new
        if(show_fixes)
            
        end
    end
    
    
    
    prev_peak_locations = peak_locations{slice_no};
end

for slice_no = 1:length(smoothed_slice_in_distr)
   
    % calculate transform parameters;

    %loc of 1.0: 
    coordinate_transform_parameters{slice_no}(1) = peak_locations{slice_no}(2);
    %loc of 2.0
    coordinate_transform_parameters{slice_no}(2) = peak_locations{slice_no}(4);
    %range between 1.0 and 2.0
    unit_transform = peak_locations{slice_no}(4) - peak_locations{slice_no}(2);

    %transform to new co-ordinate system
    origin = peak_locations{slice_no}(2) - 2*unit_transform;
    for p = 1:5
        peak_locations{slice_no}(p) = (peak_locations{slice_no}(p)-origin)/unit_transform;
    end

end


end

