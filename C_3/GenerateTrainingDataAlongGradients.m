addpath('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\PotentiallyUsefulLibraries\NIfTI_20140122');
addpath('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\PotentiallyUsefulLibraries');

%%Load DICOM (intensity data)
[Slices, ct3D, CT_dimension_spacing] = LoadDICOM(Configuration.ExampleDicomDataPath, false);


x_ct = (CT_dimension_spacing(1)/2):CT_dimension_spacing(1):CT_dimension_spacing(1)*(size(ct3D,1)-(CT_dimension_spacing(1)/2));
y_ct = (CT_dimension_spacing(2)/2):CT_dimension_spacing(2):CT_dimension_spacing(2)*(size(ct3D,2)-(CT_dimension_spacing(2)/2));
z_ct = (CT_dimension_spacing(3)/2):CT_dimension_spacing(3):CT_dimension_spacing(3)*(size(ct3D,3));
%X = {x,y,z};
X_ct = {x_ct,y_ct,z_ct};

%Load .nii bone model

display = false;


nii = load_nii(Configuration.NifitSkeletonDataPath);


bone_point_cloud_dense = nii.img;
bone_point_cloud_dense(bone_point_cloud_dense ~= 0) = 1;
%FLIP to align orientation to ct3D
bone_point_cloud_dense = flip(bone_point_cloud_dense,3);
bone_point_cloud_dense = imrotate(bone_point_cloud_dense, 90);
bone_point_cloud_dense = flip(bone_point_cloud_dense,2);


%isosurface
iso = isosurface(bone_point_cloud_dense, 0.9);
ct_coord_vertices = iso.vertices;
%figure
%patch(iso);

vertices = nan(size(ct_coord_vertices));

%transform vertices to correct co-ordinates
for dp = 1:size(ct_coord_vertices,1)
    for i = 1:3
        vertices(dp,i) = X_ct{i}(round(ct_coord_vertices(dp,i)));
    end
    
   
end

%first calculate normal directions
n = isonormals(x_ct, y_ct, z_ct, ct3D, vertices);

%now magnitudes
nm = cellfun(@norm, num2cell(n, 2));


if(bone_set)
    %remove small normals <30
   [sni] = find(nm > 50);   %THIS IS NORMAL OPERATION
    
    %%%%% TESTING / LOOKING AT bad edges
     %[sni] = find(nm < 50);

    good_set_magnitudes = nm(sni);
    good_set_positions = vertices(sni,:);
    good_set_normals = n(sni,:);
    
    
else
    good_set_magnitudes = nm;
    good_set_positions = vertices;
    good_set_normals = n;
end

%remove vertices mear the edge of the ct data;
lr = 10; %this is the radius we need to leave clear


bounds = nan(3,2);
for i = 1:3
    bounds(i,1) = lr;
    bounds(i,2) = X_ct{i}(end)-lr;
    
end


edge_point_is = [];
for dp = 1:size(good_set_positions,1)
    if(    good_set_positions(dp,1) < bounds(1,1) || good_set_positions(dp,1) > bounds(1,2)...
        || good_set_positions(dp,2) < bounds(2,1) || good_set_positions(dp,2) > bounds(2,2)...
        || good_set_positions(dp,3) < bounds(3,1) || good_set_positions(dp,3) > bounds(3,2))
        edge_point_is = [edge_point_is dp];
    end
end


%remove edge points
good_set_magnitudes(edge_point_is) = [];
good_set_positions(edge_point_is,:) = [];
good_set_normals(edge_point_is,:) = [];




