function [updated_bone_map,internal_map, num_dir_bone_found_in] = InternalPointsInline(input_mask,bone_map, r, directions_which_can_be_missed)
%INTERNALPOINTS Determines if points in the input mask are internal to the
%input mask : valid area mask (i.e. points to check whether to fill)
%bone_map
%simple: true means 6-directionality, false means 14 directionality
%


%pre-allocate
%bone_dist_maps;
num_dir_bone_found_in = uint8(zeros(size(input_mask))); %this counts the number of directions bone is found in


%bone_dist_maps - how far has been searched for bone along the specific
%direction from the specified point

%bone_dir_breakdowns - 


    num_dir = 14;
    bone_dist_maps = int16(ones(14,size(bone_map,1),size(bone_map,2),size(bone_map,3)));
    nearby_bone_map_by_direction = uint8(zeros(14,size(bone_map,1),size(bone_map,2),size(bone_map,3)));
    


%pad bone map with bone of border r!
padded_bone_map = ones(size(bone_map,1)+2*r,size(bone_map,2)+2*r,size(bone_map,3)+2*r);
padded_bone_map(r+1:size(bone_map,1)+r,r+1:size(bone_map,2)+r,r+1:size(bone_map,3)+r) = bone_map;
padded_bone_map=uint8(padded_bone_map);

%for each point in the imput mask
test_point_list = find(input_mask);
%remove bone points
bone_point_list = find(bone_map);
test_point_list(ismember(test_point_list,bone_point_list)) = [];
%i2s
[s1,s2,s3] = ind2sub(size(input_mask),test_point_list);


total_num_test_points = length(test_point_list);

%iterate through test points
%tic
for i_tp = 1:length(test_point_list)
%     if(mod(i_tp,round(length(test_point_list)/100))==0)
%         toc
%         fprintf('|');
%         tic
%     end
    cs1 = s1(i_tp);
    cs2 = s2(i_tp);
    cs3 = s3(i_tp);
    %cs = [cs1,cs2,cs3];
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
            nearby_bone_map_by_direction(1,cs1:cs1+cur_val-1,cs2,cs3) = 1;
            %num_bone_directions(cs1:cs1+cur_val-1,cs2,cs3) = num_bone_directions(cs1:cs1+cur_val-1,cs2,cs3)+1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
%     if(cur_val==r)
%         bone_dist_maps(1,cs1:cs1+cur_val-1,cs2,cs3) = r:-1:1;
%     end

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
            nearby_bone_map_by_direction(2,cs1-(cur_val-1):cs1,cs2,cs3) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
%     if(cur_val==r)
%         bone_dist_maps(2,cs1-(cur_val-1):cs1,cs2,cs3) = r:-1:1;
%     end
    
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
            nearby_bone_map_by_direction(3,cs1,cs2:cs2+(cur_val-1),cs3) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
%     if(cur_val==r)
%         bone_dist_maps(3,cs1,cs2:cs2+(cur_val-1),cs3) = r:-1:1;
%     end
    
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
            nearby_bone_map_by_direction(4,cs1,cs2-(cur_val-1):cs2,cs3) =1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
%     if(cur_val==r)
%         bone_dist_maps(4,cs1,cs2-(cur_val-1):cs2,cs3) = r:-1:1;
%     end
    
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
            nearby_bone_map_by_direction(5,cs1,cs2,cs3:cs3+(cur_val-1)) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
%     if(cur_val==r)
%         bone_dist_maps(5,cs1,cs2,cs3:cs3+(cur_val-1)) = r:-1:1;
%     end
    
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
            nearby_bone_map_by_direction(6,cs1,cs2,cs3-(cur_val-1):cs3) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
%     if(cur_val==r)
%         bone_dist_maps(6,cs1,cs2,cs3-(cur_val-1):cs3) = r:-1:1;
%     end


      %UP
      %NE
      cur_val = bone_dist_maps(7,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1+cur_val,pcs2+cur_val,pcs3+cur_val))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            
            update_points = sub2ind(size(bone_dist_maps),7*ones(cur_val,1)',...
                                    cs1:cs1+(cur_val-1),cs2:cs2+(cur_val-1),cs3:cs3+(cur_val-1));
            
            
            bone_dist_maps(update_points) = r+1;
            nearby_bone_map_by_direction(update_points) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    %NW
      cur_val = bone_dist_maps(8,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1+cur_val,pcs2+cur_val,pcs3-cur_val))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            
            update_points = sub2ind(size(bone_dist_maps),8*ones(cur_val,1)',...
                                    cs1:cs1+(cur_val-1),cs2:cs2+(cur_val-1),cs3:-1:cs3-(cur_val-1));
            
            
            bone_dist_maps(update_points) = r+1;
            nearby_bone_map_by_direction(update_points) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    %SE
      cur_val = bone_dist_maps(9,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1+cur_val,pcs2-cur_val,pcs3+cur_val))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            
            update_points = sub2ind(size(bone_dist_maps),9*ones(cur_val,1)',...
                                    cs1:cs1+(cur_val-1),cs2:-1:cs2-(cur_val-1),cs3:cs3+(cur_val-1));
            
            
            bone_dist_maps(update_points) = r+1;
            nearby_bone_map_by_direction(update_points) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    %SW
      cur_val = bone_dist_maps(10,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1+cur_val,pcs2-cur_val,pcs3-cur_val))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            
            update_points = sub2ind(size(bone_dist_maps),10*ones(cur_val,1)',...
                                    cs1:cs1+(cur_val-1),cs2:-1:cs2-(cur_val-1),cs3:-1:cs3-(cur_val-1));
            
            
            bone_dist_maps(update_points) = r+1;
            nearby_bone_map_by_direction(update_points) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    
      
      %DOWN
      
      %NE
      cur_val = bone_dist_maps(11,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1-cur_val,pcs2+cur_val,pcs3+cur_val))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            
            update_points = sub2ind(size(bone_dist_maps),11*ones(cur_val,1)',...
                                    cs1:-1:cs1-(cur_val-1),cs2:cs2+(cur_val-1),cs3:cs3+(cur_val-1));
            
            
            bone_dist_maps(update_points) = r+1;
            nearby_bone_map_by_direction(update_points) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    %NW
      cur_val = bone_dist_maps(12,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1-cur_val,pcs2+cur_val,pcs3-cur_val))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            
            update_points = sub2ind(size(bone_dist_maps),12*ones(cur_val,1)',...
                                    cs1:-1:cs1-(cur_val-1),cs2:cs2+(cur_val-1),cs3:-1:cs3-(cur_val-1));
            
            
            bone_dist_maps(update_points) = r+1;
            nearby_bone_map_by_direction(update_points) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    %SE
      cur_val = bone_dist_maps(13,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1-cur_val,pcs2-cur_val,pcs3+cur_val))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            
            update_points = sub2ind(size(bone_dist_maps),13*ones(cur_val,1)',...
                                    cs1:-1:cs1-(cur_val-1),cs2:-1:cs2-(cur_val-1),cs3:cs3+(cur_val-1));
            
            
            bone_dist_maps(update_points) = r+1;
            nearby_bone_map_by_direction(update_points) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    %SW
      cur_val = bone_dist_maps(14,cs1,cs2,cs3);
    while(cur_val<r)
        %test for bone
        if(padded_bone_map(pcs1-cur_val,pcs2-cur_val,pcs3-cur_val))
            %if bone set val and break
            %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
            %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
            %update points along line
            
            update_points = sub2ind(size(bone_dist_maps),14*ones(cur_val,1)',...
                                    cs1:-1:cs1-(cur_val-1),cs2:-1:cs2-(cur_val-1),cs3:-1:cs3-(cur_val-1));
            
            
            bone_dist_maps(update_points) = r+1;
            nearby_bone_map_by_direction(update_points) = 1;
            break;
        else
            %if not bone increment
            cur_val = cur_val+1;
        end
    end
    
    
    
end


for i_d = 1:num_dir
    num_dir_bone_found_in = num_dir_bone_found_in+ squeeze(nearby_bone_map_by_direction(i_d,:,:,:));
end
%imshow3D(num_bone_directions,[0 max(max(max(num_bone_directions)))]);
%imshow3D(squeeze(bone_dir_breakdowns(1,:,:,:)),[0 1]);
internal_map = zeros(size(num_dir_bone_found_in));
internal_map(num_dir_bone_found_in>=(num_dir-directions_which_can_be_missed)) = 1;
internal_map(input_mask==0) = 0;
updated_bone_map = bone_map;
updated_bone_map(internal_map==1) = 1;


end

