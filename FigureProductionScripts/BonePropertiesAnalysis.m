%%% Script to analyse the properties of bone seen in a CT scan (G_1)

addpath('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\PotentiallyUsefulLibraries\NIfTI_20140122');

%%Load DICOM (intensity data)
[Slices, ct3D, CT_dimension_spacing] = LoadDICOM(Configuration.ExampleDicomDataPath, false);

%% Load .nii bone model
nii = load_nii(Configuration.NifitSkeletonDataPath);
%view_nii(nii);

bone_point_cloud_dense = nii.img;
bone_point_cloud_dense(bone_point_cloud_dense ~= 0) = 1;
%FLIP to align orientation to ct3D
bone_point_cloud_dense = flip(bone_point_cloud_dense,3);
bone_point_cloud_dense = imrotate(bone_point_cloud_dense, 90);
bone_point_cloud_dense = flip(bone_point_cloud_dense,2);

sliceI = 2;
%check
%figure
%imshow(bone_point_cloud_dense(:,:,sliceI), [0 1]);
figure
slice = ct3D(:,:,sliceI);
imshow(slice, [0 max(max(slice))]);
hold on
highlight = bone_point_cloud_dense(:,:,sliceI).*max(max(slice));
highlight(highlight == 0) = nan;
transparency_map = ones(size(highlight));
transparency_map(0 == (highlight)) = 0;

g = imshow(highlight, [0 max(max(slice))]);
set(g, 'AlphaData', transparency_map);
%% properties we can anaylse

%intensity distribution
%apply mask
masked_intensity_map = nan(size(ct3D));
masked_intensity_map(bone_point_cloud_dense == 1) = ct3D(bone_point_cloud_dense == 1);
%THIS DRAWS THE INTENSITY HISTOGRAM FOR THE BONE
figure
hist_bone = histogram(masked_intensity_map);


non_bone_masked_intensity_map = nan(size(ct3D));
non_bone_masked_intensity_map(bone_point_cloud_dense == 0) = ct3D(bone_point_cloud_dense == 0);
figure
hist_non_bone = histogram(non_bone_masked_intensity_map);


%volume (by slice -> position along body)

%first calc by slice
zero_masked_bone_intensity_map = masked_intensity_map;
zero_masked_bone_intensity_map(isnan(zero_masked_bone_intensity_map)) = 0;

for si = 1:size(ct3D,3)
    bone_vol(si) = sum(sum(zero_masked_bone_intensity_map(:,:,si)));
end

%brief look
figure
plot(bone_vol)

%then average across nearby slices (moving average)
bone_vol_means = movmean(bone_vol, 20);


%then plot
figure
plot(bone_vol_means);


%edges / gradients
%get edge map
%blur then isosurface at 0.5
smoothed = imgaussfilt3(zero_masked_bone_intensity_map);
figure
slice = smoothed(:,:,sliceI);
imshow(slice, [0 max(max(slice))]);

%isosurface
iso = isosurface(bone_point_cloud_dense, 0.9);
%figure
%patch(iso);

%first calculate normal directions
n = isonormals(ct3D,iso.vertices);

%now magnitudes
nm = cellfun(@norm, num2cell(n, 2));

%quick show
figure
edges = 0:1:max(nm);
hist_surface_normal_magnitudes = histogram(nm, edges);

%% investigate small magnitude normals
small_normals = nm;
small_normals(small_normals > 20) = nan;
while(any(~isnan(nm(:))))
    [cur_normal_mag, i] = min(nm);
    min_vertex = iso.vertices(i,:);

    figure
    slice_to_show = ct3D(:,:,round(min_vertex(3)));
    imshow(slice_to_show, [0, max(max(slice_to_show))]);
    hold on;
    %box around location
    viscircles([min_vertex(1), min_vertex(2)], [5]);
    
    waitforbuttonpress;
    close
    
    %purge checked normal
    %nm(i) = nan;
end

%find small normals <30
[sni] = find(nm < 30);

%nm
%iso.vertices
%n : 3-part normals
bad_set_magnitudes = nm(sni);
bad_set_positions = iso.vertices(sni,:);
bad_set_normals = n(sni,:);

%check to see distribution of "bad" edges by 'slice'
figure
histogram(bad_set_positions(:,3));

%most occur around 200 (pelvic region)
%least in legs (which are easiest)


%remove small normals <30
[sni] = find(nm > 30);

good_set_magnitudes = nm(sni);
good_set_positions = iso.vertices(sni,:);
good_set_normals = n(sni,:);



%split into normals and magnitudes by segment along body 
%1:5:561
segment_edges = 1:5:461;
verticesy = iso.vertices(:,3);
for segment = 2:length(segment_edges)
    indices = find(verticesy > segment_edges(segment-1) & verticesy <= segment_edges(segment));
    segment_vertices{segment-1} = iso.vertices(indices,:);
    
    segment_normal_magnitudes{segment-1} = nm(indices);
end


for segment = 1:length(segment_edges)-1
    figure
    
    segment_normal_magnitude_histograms{segment} = histogram(segment_normal_magnitudes{segment});
    
    title(strcat('segment: ' ,num2str(segment_edges(segment)) , ' to ' , num2str(segment_edges(segment+1))));
   
end


%could we check the alignment of bone by assuring correct surface normal
%distribution?

%% Investigate the locality of the "good" edges set

% first make sure we have all of the "unit" information

% volume ineterpolation from old data to new data volume

gridVectCells = {0:CT_dimension_spacing(1):CT_dimension_spacing(1)*(size(ct3D,1)-1), ...
                                0:CT_dimension_spacing(2):CT_dimension_spacing(2)*(size(ct3D,2)-1), ...
                                0:CT_dimension_spacing(3):CT_dimension_spacing(3)*(size(ct3D,3)-1)};

gridIntep = griddedInterpolant(gridVectCells, ct3D, 'makima');

%dataset to test:
%good set

%parameters
r = 5; %11x11x11 with mm spacing
spacing = 1;

for dp = 1:length(good_set_positions,1)
    %interp to 2r+1 by 2r+1 by 2r+1 cube around the vertex
    unit_normal = norm(good_set_normals(dp,:));
    
end




