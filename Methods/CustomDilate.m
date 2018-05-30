function [dilated_image, new_points_image] = CustomDilate(binary_input_map, dilation_radius)
%CUSTOMDILATE Custom version of imdilate for 3D


dilated_image = binary_input_map;
bordered_dil_image = zeros(size(dilated_image,1)+2*dilation_radius,...
                            size(dilated_image,2)+2*dilation_radius,...
                            size(dilated_image,3)+2*dilation_radius);
bordered_dil_image(dilation_radius+1:dilation_radius+size(dilated_image,1),...
                    dilation_radius+1:dilation_radius+size(dilated_image,2),...
                    dilation_radius+1:dilation_radius+size(dilated_image,3)) = dilated_image;

points = find(bordered_dil_image);
[bs1,bs2,bs3] = ind2sub(size(bordered_dil_image),points);


dil_r_21 = dilation_radius*2 + 1;
dil_mask = ones(dil_r_21,dil_r_21,dil_r_21);
 

new_dil_bordered_img = zeros(size(bordered_dil_image));

num_points = length(bs1);
pre_allocate_length = num_points*26;
extended_sub_list = nan(pre_allocate_length,3);

row_pointer = 1;
tic
for d1 = -1:1:1
    for d2 = -1:1:1
        for d3 = -1:1:1
            extended_sub_list(row_pointer:row_pointer+num_points-1,1) = bs1+d1;
            extended_sub_list(row_pointer:row_pointer+num_points-1,2) = bs2+d2;
            extended_sub_list(row_pointer:row_pointer+num_points-1,3) = bs3+d3;
            row_pointer = row_pointer + num_points;
        end
    end
end
toc
tic
dilated_points = sub2ind(size(bordered_dil_image),extended_sub_list(:,1),extended_sub_list(:,2),extended_sub_list(:,3));
new_dil_bordered_img(dilated_points) = 1;
toc






% 
% for i_p = 1:length(points)
%     
%     bordered_dil_image(bs1(i_p)-dilation_radius:bs1(i_p)+dilation_radius,...
%                         bs2(i_p)-dilation_radius:bs2(i_p)+dilation_radius,...
%                         bs3(i_p)-dilation_radius:bs3(i_p)+dilation_radius) = dil_mask;
%     
% end


dilated_image = new_dil_bordered_img(dilation_radius+1:dilation_radius+size(dilated_image,1),...
                    dilation_radius+1:dilation_radius+size(dilated_image,2),...
                    dilation_radius+1:dilation_radius+size(dilated_image,3));
new_points_image = dilated_image;
new_points_image(binary_input_map==1) = 0;


end

