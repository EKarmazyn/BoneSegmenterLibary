function [valid_mask] = CreateValidAreaMask(ct3D)
%CREATEVALIDAREAMASK Summary of this function goes here

mask = ones(size(ct3D));
ct3D(400:end,:,:) = 0;


dia_mask = strel('diamond',256);
dia_mask = dia_mask.Neighborhood(1:end-1,1:end-1);

strel1 = strel('diamond',30);
dia_mask = imdilate(dia_mask,strel1);
%edge map to skin to mask
for k = 1:size(ct3D,3)
    cur_slice = ct3D(:,:,k);
    edge_map = edge(cur_slice);
    %i = find(edge_map',1);
    %[J,I] = ind2sub(size(cur_slice),i);
    CC = bwconncomp(edge_map);
    longest = 0;
    %it_longest;
    for it = 1:length(CC.PixelIdxList)
        if(length(CC.PixelIdxList{it})>longest)
            longest = length(CC.PixelIdxList{it});
            it_longest = it;
        end
    end
    lon_edge  = CC.PixelIdxList{it_longest};
    %plot longest path only
    slice_mask = zeros(size(cur_slice));
    slice_mask(lon_edge) = 1;
    
    %dilate the path
    %add line at 400;
    slice_mask(400,:) = 1;
    %add line at edges;
    slice_mask(:,1) = 1;
    slice_mask(:,end) = 1;
    se = strel('diamond',5);
    se2 = strel('diamond',8);
    
    slice_mask = imdilate(slice_mask,se);
    
    
    
    
    %fill the path
    slice_mask = imfill(slice_mask, 'holes');
    
    
    %erode
    slice_mask = imerode(slice_mask,se2);
    %slice_mask = imdilate(slice_mask,se3);
    
    
    slice_mask = slice_mask.*dia_mask;
    
    %done
    mask(:,:,k) = dia_mask;
    
end
valid_mask = uint8(mask);
end

