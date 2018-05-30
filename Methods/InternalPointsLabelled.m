function [updated_bone_map,internal_map] = InternalPointsLabelled(input_mask,labelled_map, simple, r)
%INTERNALPOINTS Determines if points in the input mask are internal to the
%bone_map, uses the labelled map to only fill within CVs
%simple: true means 6-directionality, false means 14 directionality

if(nargin<3)
    simple = false;
end
if(nargin<4)
    r = 20; %distance in voxels (diags are further
end

%pre-allocate
%bone_dist_maps;
num_bone_directions = zeros(size(input_mask));

if(simple)
    num_dir = 6;
    bone_dist_maps = ones(6,size(bone_map,1),size(bone_map,2),size(bone_map,3));
    bone_dir_breakdowns = zeros(6,size(bone_map,1),size(bone_map,2),size(bone_map,3));
else
    num_dir = 14;
    bone_dist_maps = ones(14,size(bone_map,1),size(bone_map,2),size(bone_map,3));
    bone_dir_breakdowns = zeros(14,size(bone_map,1),size(bone_map,2),size(bone_map,3));
end

%pad bone map with bone of border r!
padded_bone_map = ones(size(bone_map,1)+2*r,size(bone_map,2)+2*r,size(bone_map,3)+2*r);
padded_bone_map(r+1:size(bone_map,1)+r,r+1:size(bone_map,2)+r,r+1:size(bone_map,3)+r) = bone_map;


%for each point in the imput mask
test_point_list = find(input_mask);
[s1,s2,s3] = ind2sub(size(input_mask),test_point_list);

for i_tp = 1:length(test_point_list)
    cs1 = s1(i_tp);
    cs2 = s2(i_tp);
    cs3 = s3(i_tp);
    pcs1 = cs1+r;
    pcs2 = cs2+r;
    pcs3 = cs3+r;
    %d1+
    cur_val = bone_dist_maps(1,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1+cur_val,pcs2,pcs3))
            %if bone set val and break
            %bone_dist_maps(1,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            bone_dist_maps(1,cs1:cs1+cur_val-1,cs2,cs3) = r+1;
            bone_dir_breakdowns(1,cs1:cs1+cur_val-1,cs2,cs3) = 1;
            %num_bone_directions(cs1:cs1+cur_val-1,cs2,cs3) = num_bone_directions(cs1:cs1+cur_val-1,cs2,cs3)+1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end

    %d1-
    cur_val = bone_dist_maps(2,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1-cur_val,pcs2,pcs3))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            bone_dist_maps(2,cs1-(cur_val-1):cs1,cs2,cs3) = r+1;
            bone_dir_breakdowns(2,cs1-(cur_val-1):cs1,cs2,cs3) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    %d2+
    cur_val = bone_dist_maps(3,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1,pcs2+cur_val,pcs3))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            bone_dist_maps(3,cs1,cs2:cs2+(cur_val-1),cs3) = r+1;
            bone_dir_breakdowns(3,cs1,cs2:cs2+(cur_val-1),cs3) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end

    %d2-
    cur_val = bone_dist_maps(4,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1,pcs2-cur_val,pcs3))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            bone_dist_maps(4,cs1,cs2-(cur_val-1):cs2,cs3) = r+1;
            bone_dir_breakdowns(4,cs1,cs2-(cur_val-1):cs2,cs3) =1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    
    %d3+
    cur_val = bone_dist_maps(5,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1,pcs2,pcs3+cur_val))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            bone_dist_maps(5,cs1,cs2,cs3:cs3+(cur_val-1)) = r+1;
            bone_dir_breakdowns(5,cs1,cs2,cs3:cs3+(cur_val-1)) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    %d3-
    cur_val = bone_dist_maps(6,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1,pcs2,pcs3-cur_val))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            bone_dist_maps(6,cs1,cs2,cs3-(cur_val-1):cs3) = r+1;
            bone_dir_breakdowns(6,cs1,cs2,cs3-(cur_val-1):cs3) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    
end

for i_d = 1:num_dir
    num_bone_directions = num_bone_directions+ squeeze(bone_dir_breakdowns(i_d,:,:,:));
end
%imshow3D(num_bone_directions,[0 max(max(max(num_bone_directions)))]);
%imshow3D(squeeze(bone_dir_breakdowns(1,:,:,:)),[0 1]);
internal_map = zeros(size(num_bone_directions));
internal_map(num_bone_directions==num_dir) = 1;
updated_bone_map = bone_map;
updated_bone_map(internal_map==1) = 1;


end

