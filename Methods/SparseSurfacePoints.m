function [surface_points, surface_points_original_indexing, surface_point_normals] = SparseSurfacePoints(binary_map, axis_vectors, nominal_point_spacing)
%SPARSESURFACEPOINTS Generates a set of sparse surface points from a 3d
%binary map

if(nargin<3)
    nominal_point_spacing = 0.005; %1 point per square cm
end

%downsample
for D = 1:3
    spacing(D) = (abs(axis_vectors{D}(2) -  axis_vectors{D}(1)))/1000; %convert from mm to m
    step(D) = nominal_point_spacing/spacing(D);
    downsampled_axis{D} = round(1:step(D):size(binary_map,D));
end




downsampled_map = binary_map(downsampled_axis{1},downsampled_axis{2},downsampled_axis{3});

%new_axis_vectors
for D = 1:3
    axN{D} = axis_vectors{D}(downsampled_axis{D});
    
end



%isosurface
iso = isosurface(axN{1},axN{2},axN{3},downsampled_map,0.99);
surface_points = iso.vertices;
iso2 = isosurface(downsampled_map,0.99);
surface_points_original_indexing = round(iso2.vertices);

for D = 1:3
    surface_points_original_indexing(:,D) = downsampled_axis{D}(surface_points_original_indexing(:,D));
end


%convert downsampled surface points to more accurate positions 



%flip correct
surface_points_temp = surface_points;
surface_points(:,1) = surface_points_temp(:,2);
surface_points(:,2) = surface_points_temp(:,1);
surface_ind_points_temp = surface_points_original_indexing;
surface_points_original_indexing(:,1) = surface_ind_points_temp(:,2);
surface_points_original_indexing(:,2) = surface_ind_points_temp(:,1);
end

