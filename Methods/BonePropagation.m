function [bone_intensity_map] = BonePropagation(intensity_map,bone_binary_map, minimum_bone_intensity)
%BONEPROPAGATION 

bone_intensity_map = zeros(size(intensity_map));
bone_intensity_map(bone_binary_map == true) = intensity_map(bone_binary_map == true);

%apply a filter to calculate average local bone intensity
% r = 6;
% b = r+1;
% SE = strel('octagon',r);
% SE_ROI = strel('square',3);
% bone_locality_kernel = double(SE.Neighborhood);
% bone_locality_kernel(r+1,r+1) = 0;
% %nearby_bone_intensity_map = imfilter(bone_intensity_map,bone_locality_kernel);

%nearby_bone_intensity_map = zeros(size(bone_intensity_map));

%bordered_bone_intensity_map = zeros(size(intensity_map,1)+2*(b-1),size(intensity_map,1)+2*(b-1));
%bordered_bone_intensity_map((b):(end-(r)),(b):(end-(r))) = bone_intensity_map;

%bordered_bone_binary_map = zeros(size(intensity_map,1)+2*(b-1),size(intensity_map,1)+2*(b-1));
%bordered_bone_binary_map((b):(end-(r)),(b):(end-(r))) = bone_binary_map;

%still_change = true;
%iters = 0;

%nearby_bone_intensity_map_wo_divide = imfilter(bone_intensity_map,bone_locality_kernel);
%nearby_bone_intensity_map_num = imfilter(bone_binary_map,bone_locality_kernel);
%nearby_bone_intensity_map = nearby_bone_intensity_map_wo_divide./nearby_bone_intensity_map_num;



potential_bone_map = zeros(size(intensity_map));
potential_bone_map(intensity_map > minimum_bone_intensity) = 1;

CC = bwconncomp(potential_bone_map);

%calculate gradient of the potential bone map
[FX,FY] = gradient(intensity_map);
GRAD = ((FX.^2).*(FY.^2)).^(1/2);

edge_map = edge(intensity_map, 'Sobel', 80);

figure
imshow(edge_map);

smoothed_intensity_map = imgaussfilt(intensity_map,4);
figure
imshow(smoothed_intensity_map,[0 max(max(smoothed_intensity_map))]);


edge_map = edge(smoothed_intensity_map, 'Sobel');

figure
imshow(edge_map);
%remove areas with too large of a gradient (







% 
% 
% still_change = true;
% figure
% while(still_change)
%     imshow(bone_binary_map);
%     iters = iters + 1;
%     
%     %generate ROI mask
%     dil_bone_binary_map = imdilate(bone_binary_map,SE_ROI);
%     roi_mask = xor(dil_bone_binary_map, bone_binary_map);
%     
%     
%     
%     nearby_bone_intensity_map_wo_divide = roifilt2(bone_locality_kernel, bone_intensity_map, roi_mask);
%     nearby_bone_intensity_map_num = roifilt2(bone_locality_kernel,bone_binary_map,roi_mask);
%     nearby_bone_intensity_map(roi_mask) = nearby_bone_intensity_map_wo_divide(roi_mask)./nearby_bone_intensity_map_num(roi_mask);
%     
%     %%%
%     a = nearby_bone_intensity_map(roi_mask)*Configuration.BonePropagation_NearbyBoneFactor;
%     
%     bone_binary_map(roi_mask) = (intensity_map(roi_mask)>minimum_bone_intensity  & a < intensity_map(roi_mask));
%     bone_intensity_map(bone_binary_map == 1) = intensity_map(bone_binary_map == 1);
%     
% end
% 
% figure
% 
% while(still_change)
%     %close
%     
%     imshow(bordered_bone_intensity_map);
%     
%     iters = iters+1;
%     still_change = false;
%     for i = 1:size(nearby_bone_intensity_map,1)
%         bi = i+r;
%         for j = 1:size(nearby_bone_intensity_map,2)
%             bj = j+r;
%             %for k = 1:size(nearby_bone_intensity_map,3)
%                 %bk = b+k;
%                 nearby_count = sum(sum(bordered_bone_binary_map(bi-r:bi+r, bj-r:bj+r)));
%                 if(nearby_count ~= 0)
%                     nearby = bordered_bone_intensity_map(bi-r:bi+r, bj-r:bj+r);
% 
%                     val = sum(sum(nearby))/nearby_count;
%                     nearby_bone_intensity_map(i,j) = val;
%                     if(intensity_map(i,j) > Configuration.BonePropagation_NearbyBoneFactor*val)
%                         bordered_bone_binary_map(bi,bj) = 1;
%                         bordered_bone_intensity_map(bi,bj) = intensity_map(i,j);
%                         still_change = true;
%                     end
%                 end
% 
%             %end
%         end
%     end
%     
%     title(num2str(iters));
% end

end

