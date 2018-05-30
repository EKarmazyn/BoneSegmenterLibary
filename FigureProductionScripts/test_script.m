[Slices, ct3D] = LoadDICOM(Configuration.Example2DicomDataPath, false);


%%
binary_map = zeros(size(ct3D));
V = 500;
binary_map(ct3D > V) = 1;

%%

si = 1;
slice = binary_map(:,:,si);
%CC = bwconncomp();

BW = edge(binary_map(:,:,si));
figure
imshow(BW);


SE = strel('octagon', 3);
BW2 = imdilate(BW, SE);
figure
imshow(BW2);

BW3 = imerode(BW2, SE);
figure
imshow(BW3);




[B,L] = bwboundaries(BW3, 'noholes');
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
end

longest_boundary_index = [];
longest_boundary_length = 0;
for k = 1:length(B)
    boundary_length = length(B{k});
    if(boundary_length > longest_boundary_length)
        longest_boundary_length = boundary_length;
        longest_boundary_index = k;
    end
end
body_boundary = B{longest_boundary_index};
plot(body_boundary(:,2), body_boundary(:,1), 'g', 'LineWidth', 2)

bin_slice = zeros(size(slice));
bin_slice(L == longest_boundary_index) = 1;

bin_slice = imfill(bin_slice, 'holes');

% flesh_binary = zeros(size(slice));
% flesh_binary(ct3D(:,:,si) < 1000) = 1;
% 
% 
% 
% flesh_binary(bin_slice == 0) = 0;
% 
% figure
% imshow(flesh_binary);

%%
intensity_map = ct3D(:,:,si);
valid_areas_map = bin_slice;
[bone_seed_map, minimum_bone_intensity] = FindBoneSeeds(intensity_map, valid_areas_map);


bone_probability_map = nan(size(intensity_map));
bone_probability_map(bone_seed_map == 1) = 1;
fixed_map = false(size(intensity_map));

[bone_intensity_map] = BonePropagation(intensity_map,bone_seed_map, minimum_bone_intensity);