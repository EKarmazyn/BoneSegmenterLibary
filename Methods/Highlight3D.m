function Highlight3D(fuzzy_bone_map, connected_volume_array_red, connected_volume_array_green, connected_volume_array_blue)
%HIGHLIGHT3D uses imshow3d and highlights the passed connected volume

 
color_cc_map = zeros(size(fuzzy_bone_map,1),size(fuzzy_bone_map,2),size(fuzzy_bone_map,3),3);


    
    for cv = 1:length(connected_volume_array_red)
        connected_volume = connected_volume_array_red(cv);
        all_points_ind_list = connected_volume.IndVoxelList;

        %ind 2 sub
        [co_sub{1},co_sub{2},co_sub{3}] = ind2sub(size(fuzzy_bone_map),all_points_ind_list);


        esub_add{4} = zeros(length(co_sub{1}),1) + 1; %+3blue, +1 is red, +2 green

        %sub to ind

        co_ind_add = sub2ind(size(color_cc_map),co_sub{1},co_sub{2},co_sub{3},esub_add{4});

        %update color map

        color_cc_map(co_ind_add) = 1;
    end
    
    for cv = 1:length(connected_volume_array_green)
        connected_volume = connected_volume_array_green(cv);
        all_points_ind_list = connected_volume.IndVoxelList;

        %ind 2 sub
        [co_sub{1},co_sub{2},co_sub{3}] = ind2sub(size(fuzzy_bone_map),all_points_ind_list);


        esub_add{4} = zeros(length(co_sub{1}),1) + 2; %+3blue, +1 is red, +2 green

        %sub to ind

        co_ind_add = sub2ind(size(color_cc_map),co_sub{1},co_sub{2},co_sub{3},esub_add{4});

        %update color map

        color_cc_map(co_ind_add) = 1;
    end
    
    for cv = 1:length(connected_volume_array_blue)
        connected_volume = connected_volume_array_blue(cv);
        all_points_ind_list = connected_volume.IndVoxelList;

        %ind 2 sub
        [co_sub{1},co_sub{2},co_sub{3}] = ind2sub(size(fuzzy_bone_map),all_points_ind_list);


        esub_add{4} = zeros(length(co_sub{1}),1) + 3; %+3blue, +1 is red, +2 green

        %sub to ind

        co_ind_add = sub2ind(size(color_cc_map),co_sub{1},co_sub{2},co_sub{3},esub_add{4});

        %update color map

        color_cc_map(co_ind_add) = 1;
    end
    
    imshow3D(color_cc_map);

end

