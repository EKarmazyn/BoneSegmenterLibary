function [surface_points, surface_points_original_indexing, surface_point_normals] = SparseSurfacePoints2(binary_map, axis_vectors, nominal_point_spacing, edge_margin)
%SPARSESURFACEPOINTS2 Different implementation of sparse surface point
%finding. Looks slice by slice. Finds 
DEBUG = false;
if(nargin<4)
    edge_margin = 0.02; %20mm
end
if(nargin<3)
    nominal_point_spacing = 0.005; %1 point per square cm
end
%calculate how many points we need
%first volume
%downsample
for D = 1:3
    spacing(D) = (abs(axis_vectors{D}(2) -  axis_vectors{D}(1)))/1000; %convert from mm to m
end

%calc margins
for D = 1:3
    ind_margin(D) = ceil(edge_margin/spacing(D));
    ind_margin_min(D) = 1+ind_margin(D);
    ind_margin_max(D) = length(axis_vectors{D})-ind_margin(D);
end


vol = sum(sum(sum(binary_map)))*(spacing(1)*spacing(2)*spacing(3));
%find r:
r = ((3*vol)/(4*pi))^(1/3);
%as sphere surface area
surf_area = vol*3/r;
num_points = surf_area/(nominal_point_spacing^2);

num_slices = size(binary_map,3);


points_per_slice = num_points/num_slices;


surface_points = []; %x,y,z
surface_points_original_indexing = [];

for i_slice = 1:size(binary_map,3)
    cur_slice = squeeze(binary_map(:,:,i_slice));
    
    
    [B,L,N,A] = bwboundaries(cur_slice, 'noholes');
    all_points = cell2mat(B);
    
    
    %filter out points near margin
    valid_map = ones(size(all_points,1),1);
    for D = 1:2
        valid_map = valid_map.*((all_points(:,D)> ind_margin_min(D)) .* (all_points(:,D)< ind_margin_max(D)));
    end
    valid_map = valid_map.*((i_slice> ind_margin_min(3)) .*( i_slice< ind_margin_max(3)));
    if(~any(valid_map))
        continue;
    end
    points_left = all_points((valid_map.*(1:size(all_points,1))'),:);
    
    
    
    
    %if num points left < points_per_slice
    if(size(points_left,1) < points_per_slice)
        surface_points_original_indexing = [surface_points_original_indexing; points_left,zeros(size(points_left,1),1)+i_slice];
    else
        
        %select N randomly
        sample = round(rand(ceil(points_per_slice),1) * (size(points_left,1)-1))+1;
        sel_points = points_left(sample,:);

        surface_points_original_indexing = [surface_points_original_indexing; sel_points, zeros(size(sel_points,1),1)+i_slice];
    
        
    end
    
    
    if(DEBUG)
        
        imshow(1-cur_slice);
        hold on;
        for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
        end
        
        scatter(sel_points(:,2), sel_points(:,1),'g');
        
    end
end

%convert to real positions
surface_points = nan(size(surface_points_original_indexing));
for D = 1:3
    surface_points(:,D) = axis_vectors{D}(surface_points_original_indexing(:,D));
    
end



end

