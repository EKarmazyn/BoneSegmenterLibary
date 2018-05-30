function [peak_locations, soft_tissue_std] = PrimitiveFindPeaksInDistr(ct_distribution, display, description)
%PRIMITIVEFINDPEAKSINDISTR Summary of this function goes here


to_display=false;
if(nargin > 1)
    to_display = display;
    if(nargin > 2)
        
    else
        description = '';
    end
end

peak_locations = nan(5,1);
    

    
[peaks, locs, w, p] = findpeaks(ct_distribution);

[peaks_sorted, sort_vector] = sort(peaks, 'descend');
locs_sorted = locs(sort_vector);
w_sorted = w(sort_vector);
p_sorted = p(sort_vector);


largest_peak_prominence = p_sorted(1);
temp = p_sorted;
temp(p_sorted < Configuration.PeakAlignment_ValidPeakParameter * largest_peak_prominence) = 0;
valid_peaks = find(temp);

num_valid_peaks = length(valid_peaks); %hopefully 5!
new_valid_peaks = [];


%check peaks arent on top of each other

if(num_valid_peaks > 5)
    % crop to top 5 peaks (by prominence value) OLD
    % NEED TO CROP TO 2 in region A and 3 in region B --------- HERE 

    %first find location of last valid peak
    %split range into thirds, two max in first third, 3 max in second
    loc_furthest_valid_peak = max(locs_sorted(valid_peaks));
    early_peaks = [];
    late_peaks = [];
    for vp = 1:num_valid_peaks
        peak_loc = locs_sorted(valid_peaks(vp));
       
        if(peak_loc > loc_furthest_valid_peak*(1/2))
            late_peaks = [late_peaks valid_peaks(vp)];
        elseif(peak_loc < loc_furthest_valid_peak*(1/3))
            early_peaks = [early_peaks valid_peaks(vp)];
        end

    end

    if(length(early_peaks) > 2)
        prom_sorted = sort(p_sorted(early_peaks) ,'descend');
        prom_cutoff = prom_sorted(2);
        to_remove = [];
        for p_num = 1:length(early_peaks)
            if(p_sorted(early_peaks(p_num)) < prom_cutoff)
                to_remove = [to_remove early_peaks(p_num)];
            end
        end
        for rm = 1:length(to_remove)
            early_peaks(early_peaks == to_remove(rm)) = [];
        end

    end

    if(length(late_peaks) > 3)
        prom_sorted = sort(p_sorted(late_peaks) ,'descend');
        prom_cutoff = prom_sorted(3);
        to_remove = [];
        for p_num = 1:length(late_peaks)
            if(p_sorted(late_peaks(p_num)) < prom_cutoff)
                to_remove = [to_remove late_peaks(p_num)];
            end
        end
        for rm = 1:length(to_remove)
            late_peaks(late_peaks == to_remove(rm)) = [];
        end
    end

    valid_peaks = [early_peaks late_peaks];
    num_valid_peaks = length(valid_peaks);
    

else
    %5 peaks
    
end

if(num_valid_peaks < 5)

warning('Algorithm found less than 5 valid peaks in slice, likely that p1 and p2 are merged, or p3/p5 are unclear');
end





first_peak_loc = min(locs_sorted(valid_peaks));
last_peak_loc = max(locs_sorted(valid_peaks));

%two in first 1/3, 3 in last 1/2
ov_range = last_peak_loc - first_peak_loc;
regionA = [(first_peak_loc-1) (first_peak_loc + (ov_range/3))];
regionB = [(last_peak_loc-(ov_range/2)) (last_peak_loc+1)];

peaks_in_A = [];
peaks_in_B = [];
for v_peak = 1:num_valid_peaks
    peak_loc = locs_sorted(valid_peaks(v_peak));

    if(peak_loc > regionA(1) && peak_loc < regionA(2))
       peaks_in_A = [peaks_in_A valid_peaks(v_peak)]; 


    elseif(peak_loc > regionB(1) && peak_loc < regionB(2))
        peaks_in_B = [peaks_in_B valid_peaks(v_peak)]; 
    end
end

valid_locs_sorted = locs_sorted(valid_peaks);
[~, peak_loc_sort] = sort(locs_sorted(valid_peaks));



if (length(peaks_in_A) == 2 && length(peaks_in_B) == 3)
    % GOOD!
    peak_locations(1) = valid_locs_sorted(peak_loc_sort(1));
    peak_locations(2) = valid_locs_sorted(peak_loc_sort(2));
    peak_locations(3) = valid_locs_sorted(peak_loc_sort(3));
    peak_locations(4) = valid_locs_sorted(peak_loc_sort(4));
    peak_locations(5) = valid_locs_sorted(peak_loc_sort(5));





else
    % region A

    if(length(peaks_in_A) > 2)
        %purge
        [min_promin, ip] = min(p_sorted(peaks_in_A));
        peaks_in_A(ip) = [];
    end

    if (length(peaks_in_A) == 1)
        % merged
        peak_locations(1) = NaN;
        peak_locations(2) = locs_sorted(peaks_in_A(1));
    elseif (length(peaks_in_A) == 2)
        % good
        [~,loc_order]  = sort(locs_sorted(peaks_in_A));
        peak_locations(1) = locs_sorted(peaks_in_A(loc_order(1)));
        peak_locations(2) = locs_sorted(peaks_in_A(loc_order(2)));

    end

    [~,loc_order]  = sort(locs_sorted(peaks_in_B));
    %region B
    if (length(peaks_in_B) == 3)
        %peaks_in_B contains the 3 peaks, they index into the
        %initial_sorted arrays. We need to sort them into loc order

        peak_locations(3) = locs_sorted(peaks_in_B(loc_order(1)));
        peak_locations(4) = locs_sorted(peaks_in_B(loc_order(2)));
        peak_locations(5) = locs_sorted(peaks_in_B(loc_order(3)));
    elseif (length(peaks_in_B) == 2)
        %get widths of peaks in b, take normal distributions away to reveal
        %third peak.
        b1 = min(valid_locs_sorted(peaks_in_B));
        b2 = max(valid_locs_sorted(peaks_in_B));
        w1 = w(find(locs==b1));
        h1 = peaks(find(locs==b1));
        
        loc_left = b1-round(w1/2); %%maybe not halfed?
        h_left = ct_distribution(loc_left);
        frac_left = h_left/h1;
        sigmasq_1 = -(round(w1/2)^2)/(2*log(frac_left));
        w2 = w(find(locs==b2));
        h2 = peaks(find(locs==b2));
        
        loc_left = b2-round(w2/2);
        h_left = ct_distribution(loc_left);
        frac_left = h_left/h2;
        sigmasq_2 = -(round(w2/2)^2)/(2*log(frac_left));
        
        %region b
        region_start = b1-100;
        region_end = min(b2+150,length(ct_distribution));
        region = region_start:region_end;
        b_region = ct_distribution(region_start:region_end);
         
        %debug%
%         figure
%         hold on
%         plot(b_region,'r');
        %debug%
        
       
        %1
        normal_1 = @(x) h1*exp(-((x-b1)^2)/(2*sigmasq_1)) ;
        for in_region_position = 1:length(region)
            b_region(in_region_position) = b_region(in_region_position)-normal_1(region_start+in_region_position-1);
        end
        
        %debug%
        %plot(b_region,'b');
        %debug%
        
        %2
        normal_2 = @(x) h2*exp(-((x-b2)^2)/(2*sigmasq_2)) ;
        for in_region_position = 1:length(region)
            b_region(in_region_position) = b_region(in_region_position)-normal_2(region_start+in_region_position-1);
        end
        
        %debug%
        %plot(b_region,'g');
        %debug%
        
        %find two largest peaks, one on the right is the 5th peak
        [new_peaks, new_locs, ~, new_p] = findpeaks(b_region);
        [sorted_p, sI] = sort(new_p, 'descend');
        best_two = new_locs(sI(1:2));
        peak_5_loc = max(best_two) + region_start;
        
        
        peak_locations(3) = b1;
        peak_locations(4) = b2;
        peak_locations(5) = peak_5_loc;
        
        
        
        %%%OLD
%         % need to find out which two peaks we have.
%         % if both are prominent it is 4 and 5,
%         % otherwise label the largest as 4
%         [larger_prominence,Il] = max(p_sorted(peaks_in_B));
%         [smaller_prominence,Is] = min(p_sorted(peaks_in_B));
%         if(smaller_prominence > Configuration.PeakAlignment_RelativeProminanceParameter*larger_prominence)
%             %4 and 5
%             peak_locations(4) = locs_sorted(peaks_in_B(loc_order(1)));
%             peak_locations(5) = locs_sorted(peaks_in_B(loc_order(2)));
%         else
%             %peak_locations{slice_no}(3) = locs_sorted(peaks_in_B(loc_order(1)));
%             peak_locations(4) = locs_sorted(peaks_in_B(Il));
%         end

    else
        %error
        %display('only one peak found in region B');
        %we have found peak 4
         peak_locations(4) = locs_sorted(peaks_in_B(1));
    end


end




if(to_display)
    %plot found peaks
    h = figure;
    hold on
    plot(ct_distribution);
    %blobs on each peak
    for peak_num = 1:5
        if(~isnan(peak_locations(peak_num)))
            %viscircles([peak_locations{slice_no}(peak_num) smoothed_slice_in_distr{slice_no}(peak_locations{slice_no}(peak_num))], 5); 
            text(peak_locations(peak_num), ct_distribution(peak_locations(peak_num)), strcat('p',num2str(peak_num)));
        end
    end
    title(description);

    waitfor(h);

end

%get std for right soft tissue peak (4)
%peak_locations(4)
soft_tissue_id = find(locs==peak_locations(4));
width = w(soft_tissue_id);
loc_left = peak_locations(4)-round(width/2);
h_left = ct_distribution(loc_left);
frac_left = h_left/peaks(soft_tissue_id);
sigmasq = -(round(width/2)^2)/(2*log(frac_left));

soft_tissue_std = sqrt(sigmasq);
        

end

