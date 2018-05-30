function [voided, voided_eroded, voided_eroded2] = HardSeperateVolumes(ct3D, binary_map_filled, CT_dimension_spacing,xy_smoothed_thresh, z_thresh, SCALE )
%HARDSEPERATEVOLUMES aggressive seperation at end
    

    %xy_smoothed_thresh = 30;

    rows = size(binary_map_filled,1);
    cols = size(binary_map_filled,2);
    slices = size(binary_map_filled,3);
    
    %xy_thresh = 0;
    %z_thresh = 10;
    %dog_gf1 = 0.0104* SCALE.xy * 1.2611;
    %dog_gf2 = 1.9997* SCALE.xy * 1.2611;
    
    bt = 0.3;
    num_1 = round(2* (0.792968750000000/SCALE.xy));
    num_2 = round(4* (0.792968750000000/SCALE.xy));
    span = num_1;
    sps = num_2;
    h = gaussdesign(bt,span,sps);
    
    span = num_2;
    sps = num_1;
    h2 = gaussdesign(bt,span,sps);
    
    h_dog = h2-h;
    filt_length_1d = length(h_dog);
    
    x_edge_map = zeros(size(ct3D));
    offset_raw = int32(ct3D) - (1000 + (SCALE.initial_t-1175));
    
    for i_y = 1:cols
        %tic
        for i_z = 1:slices
            x_array = offset_raw(:,i_y,i_z);
            edge_v = conv(x_array,h_dog,'same');
            x_edge_map(:,i_y,i_z) = edge_v(:);
        end
        %toc
    end
    
    y_edge_map = zeros(size(ct3D));
    %offset_raw = int32(ct3D) - 1000;
    
    for i_x = 1:rows
        %tic
        for i_z = 1:slices
            x_array = squeeze(offset_raw(i_x,:,i_z));
            edge_v = conv(x_array,h_dog,'same');
            y_edge_map(i_x,:,i_z) = edge_v(:);
        end
        %toc
    end
    
    xy_sum_edge_map = (x_edge_map.^2 + y_edge_map.^2).^(1/2);
%     
%     xy_dog_sum_edge_map = zeros(size(ct3D));
%     for i_z = 1:361
%         xy_dog_sum_edge_map(:,:,i_z) = imgaussfilt(xy_sum_edge_map(:,:,i_z),dog_gf1) - ...
%             imgaussfilt(xy_sum_edge_map(:,:,i_z),dog_gf2);
%      end
    
    
    z_rel_scale = CT_dimension_spacing(3)/CT_dimension_spacing(1);
    
    
    num_1 = round(1* (1.5/SCALE.z));
    num_2 = round(2* (1.5/SCALE.z));
    
    
    bt = 0.3;
    span = num_1;
    sps = num_2;
    hz = gaussdesign(bt,span,sps);
    
    span = num_2;
    sps = num_1;
    h2z = gaussdesign(bt,span,sps);
    hz_dog = h2z-hz;
    %fvtool(hz_dog,'impulse')
    
    
    z_edge_map = zeros(size(ct3D));
    for i_x = 1:rows
        %tic
        for i_y = 1:cols
            z_array = squeeze(offset_raw(i_x,i_y,:));
            edge_v = conv(z_array,hz_dog,'same');
            z_edge_map(i_x,i_y,:) = edge_v(:);
        end
        %toc
    end
    
%     figure('Name','z_edge_map')
%     imshow3D(z_edge_map, [-500,500]);
%     
%     figure('Name','xy_dog_sum_edge_map')
%     imshow3D(xy_dog_sum_edge_map, [-500,500]);
    
    %z_edge_map >10?
    %xy_dog_sum_edge_map >0?
    
    voided = binary_map_filled;
    %probable x-y edge
    
    smoothed = imgaussfilt3(double(xy_sum_edge_map),1);
    
    
    
    voided(smoothed>=xy_smoothed_thresh) = 0;
    
    %probable z edge
    
    %smoothedz = imgaussfilt3(double(z_edge_map),1);
    SE = strel('sphere',2);
    %SE1 = strel('sphere',1);
    %z_thresh 25 is good when dilating the edges
    z_edges = z_edge_map>=z_thresh;
    %only dilate vertically! first
     z_edges_dil = imdilate(z_edges,[0 0 0;0 0 0;1 1 1]);
    
   
    %imshow3D(z_edges_dil)
    voided(z_edges_dil==1) = 0;
    %old_voided = voided;
    
    
    
    voided_eroded = imerode(voided,SE);
    voided_eroded2 = imerode(voided_eroded,SE);
    
    comb = voided;
    comb(voided_eroded==1) = 2;
    comb(voided_eroded2==1) = 3;
    
    
    
    
    
    
    %combined
%     xy_semi = zeros(size(ct3D));
%     z_semi = zeros(size(ct3D));
%     xy_semi(xy_dog_sum_edge_map>= (xy_thresh-10)) = 1;
%     z_semi(z_edge_map>=(z_thresh-10)) = 1;
%     comb = xy_semi .* z_semi;
%     
%     voided(comb==1) = 0;
%     
%     old_voided(comb==1) = 2;
%    figure
%    imshow3D(old_voided, [0 2]);
    
%     [connected_volumes2, labelled_binary_map2,all_points_ind_list2] = ConnectedVolumes(voided);
%     figure
%     imshow3D(labelled_binary_map2);
    
end

