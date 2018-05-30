function [updated_bone_map,internal_map] = InternalPoints(input_mask,bone_map, simple, r)
%INTERNALPOINTS Determines if points in the input mask are internal to the
%input mask : valid area mask (i.e. points to check whether to fill)
%bone_map
%simple: true means 6-directionality, false means 14 directionality
%

if(nargin<3)
    simple = false;
end
if(nargin<4)
    r = 20; %distance in voxels (diags are further
end

%pre-allocate
%bone_dist_maps;
num_dir_bone_found_in = uint8(zeros(size(input_mask))); %this counts the number of directions bone is found in


%bone_dist_maps - how far has been searched for bone along the specific
%direction from the specified point

%bone_dir_breakdowns - 

if(simple)
    num_dir = 6;
    bone_dist_maps = uint8(ones(6,size(bone_map,1),size(bone_map,2),size(bone_map,3)));
    nearby_bone_map_by_direction = uint8(zeros(6,size(bone_map,1),size(bone_map,2),size(bone_map,3)));
    direction_search = [1 0 0;... 
                        -1 0 0;... 
                        0 1 0;... 
                        0 -1 0;... 
                        0 0 1;... 
                        0 0 -1];
else
    num_dir = 14;
    bone_dist_maps = uint8(ones(14,size(bone_map,1),size(bone_map,2),size(bone_map,3)));
    nearby_bone_map_by_direction = uint8(zeros(14,size(bone_map,1),size(bone_map,2),size(bone_map,3)));
    direction_search = [1 0 0;... 
                        -1 0 0;... 
                        0 1 0;... 
                        0 -1 0;... 
                        0 0 1;... 
                        0 0 -1;... 
                        
                        %up
                        1 1 1;... 
                        1 1 -1;... 
                        1 -1 1;... 
                        1 -1 -1;... 
                        
                        %down
                        -1 1 1;... 
                        -1 1 -1;... 
                        -1 -1 1;... 
                        -1 -1 -1];
end

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
tic
for i_tp = 1:length(test_point_list)
    if(mod(i_tp,1000)==0)
        toc
        fprintf('|');
        tic
    end
    cs1 = s1(i_tp);
    cs2 = s2(i_tp);
    cs3 = s3(i_tp);
    cs = [cs1,cs2,cs3];
    pcs1 = cs1+r;
    pcs2 = cs2+r;
    pcs3 = cs3+r;
    
    %iterate through directions
    for i_dir = 1:num_dir
        %tStartDir = tic;
        cur_search_distance = int16(bone_dist_maps(i_dir,cs1,cs2,cs3));
        %iterate until too far, points which we know are too far are
        %already set to r+1, points which have been searched for i<r are
        %set to i
        while(cur_search_distance<=r)
            %tStartWhileLoop = tic;
            %if bone at search distance in search direction
            cur_search_point_padded = [pcs1+direction_search(i_dir,1)*cur_search_distance,...
                                pcs2+direction_search(i_dir,2)*cur_search_distance,...
                                pcs3+direction_search(i_dir,3)*cur_search_distance];
            cur_search_point_minus_1_NOT_padded = [cs1+direction_search(i_dir,1)*(cur_search_distance-1),...
                                cs2+direction_search(i_dir,2)*(cur_search_distance-1),...
                                cs3+direction_search(i_dir,3)*(cur_search_distance-1)];
            %tStartIndSub = tic;
            linspace_subs = CustomLinspace(cs,direction_search(i_dir,:),cur_search_point_minus_1_NOT_padded, size(bone_map));
            update_inds = sub2ind(size(bone_dist_maps),i_dir*ones(size(linspace_subs(:,1))),linspace_subs(:,1),linspace_subs(:,2),linspace_subs(:,3));
            %fprintf("time for Ind Sub Stuff: " + num2str(toc(tStartIndSub))+ "\n");
            %             for i_udpd = 1:3
%                 update_points{i_udpd} = %cs(i_udpd):direction_search(i_dir,i_udpd):cur_search_point_minus_1_NOT_padded(i_udpd);
%                 if(isempty(update_points{i_udpd}))
%                     update_points{i_udpd} = 1;
%                 end
% 
%             end
            if(padded_bone_map(cur_search_point_padded(1),cur_search_point_padded(2),cur_search_point_padded(3)))
                
               
                
                %update_points = [cs1:direction_search(i_dir,1):cur_search_point_minus_1_NOT_padded(1); cs2:direction_search(i_dir,2):cur_search_point_minus_1_NOT_padded(2);cs3:direction_search(i_dir,3):cur_search_point_minus_1_NOT_padded(3)];
                %update search distance of points from cur point to found
                %bone (minus 1) to r+2, so they dont get searched from again
                bone_dist_maps(update_inds) = r+2;
                %update the found nearby bone map
                nearby_bone_map_by_direction(update_inds) = 1;
                break;
            %else increase search distance
            else
                %if not bone increment
                cur_search_distance = cur_search_distance+1;
            end
            %fprintf("time for while loop: " + num2str(toc(tStartWhileLoop))+ "\n");

        end
        %check if the while loop ended due to search distance going above r
        %(and not finding bone)
        if(cur_search_distance==r+1)
            %in this case: we wish to label points along the line searched
            %with their current search distance (i.e. r:-1:1)
            %labelling with 1 does nothing as unlabelled = 1; 
%             cur_search_point_minus_1_NOT_padded2 = [cs1+direction_search(i_dir,1)*(r-1),...
%                                 cs2+direction_search(i_dir,2)*(r-1),...
%                                 cs3+direction_search(i_dir,3)*(r-1)];
            bone_dist_maps(update_inds) =  r:-1:1;
                    
            %bone_dist_maps(i_dir,cs1:cur_search_point_minus_1_NOT_padded(1), cs2:cur_search_point_minus_1_NOT_padded(2),cs3:cur_search_point_minus_1_NOT_padded(3)) = r:-1:1;
                                
        end
        %fprintf("time for 1 dir: " + num2str(toc(tStartDir))+ "\n");
    end
    
%     %d1+
%     %cur_val = bone_dist_maps(1,cs1,cs2,cs3);
%     cur_val = 1;
%     while(cur_val<r)
%         %test for bone
%         if(padded_bone_map(pcs1+cur_val,pcs2,pcs3))
%             %if bone set val and break
%             %bone_dist_maps(1,cs1,cs2,cs3) = cur_val;
%             %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
%             %update points along line
%             %bone_dist_maps(1,cs1:cs1+cur_val-1,cs2,cs3) = r+1;
%             nearby_bone_map_by_direction(1,cs1:cs1+cur_val-1,cs2,cs3) = 1;
%             %num_bone_directions(cs1:cs1+cur_val-1,cs2,cs3) = num_bone_directions(cs1:cs1+cur_val-1,cs2,cs3)+1;
%             break;
%         else
%             %if not bone increment
%             cur_val = cur_val+1;
%         end
%     end
% 
%     %d1-
%     cur_val = bone_dist_maps(2,cs1,cs2,cs3);
%     while(cur_val<r)
%         %test for bone
%         if(padded_bone_map(pcs1-cur_val,pcs2,pcs3))
%             %if bone set val and break
%             %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
%             %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
%             %update points along line
%             bone_dist_maps(2,cs1-(cur_val-1):cs1,cs2,cs3) = r+1;
%             nearby_bone_map_by_direction(2,cs1-(cur_val-1):cs1,cs2,cs3) = 1;
%             break;
%         else
%             %if not bone increment
%             cur_val = cur_val+1;
%         end
%     end
%     
%     %d2+
%     cur_val = bone_dist_maps(3,cs1,cs2,cs3);
%     while(cur_val<r)
%         %test for bone
%         if(padded_bone_map(pcs1,pcs2+cur_val,pcs3))
%             %if bone set val and break
%             %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
%             %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
%             %update points along line
%             bone_dist_maps(3,cs1,cs2:cs2+(cur_val-1),cs3) = r+1;
%             nearby_bone_map_by_direction(3,cs1,cs2:cs2+(cur_val-1),cs3) = 1;
%             break;
%         else
%             %if not bone increment
%             cur_val = cur_val+1;
%         end
%     end
% 
%     %d2-
%     cur_val = bone_dist_maps(4,cs1,cs2,cs3);
%     while(cur_val<r)
%         %test for bone
%         if(padded_bone_map(pcs1,pcs2-cur_val,pcs3))
%             %if bone set val and break
%             %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
%             %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
%             %update points along line
%             bone_dist_maps(4,cs1,cs2-(cur_val-1):cs2,cs3) = r+1;
%             nearby_bone_map_by_direction(4,cs1,cs2-(cur_val-1):cs2,cs3) =1;
%             break;
%         else
%             %if not bone increment
%             cur_val = cur_val+1;
%         end
%     end
%     
%     
%     %d3+
%     cur_val = bone_dist_maps(5,cs1,cs2,cs3);
%     while(cur_val<r)
%         %test for bone
%         if(padded_bone_map(pcs1,pcs2,pcs3+cur_val))
%             %if bone set val and break
%             %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
%             %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
%             %update points along line
%             bone_dist_maps(5,cs1,cs2,cs3:cs3+(cur_val-1)) = r+1;
%             nearby_bone_map_by_direction(5,cs1,cs2,cs3:cs3+(cur_val-1)) = 1;
%             break;
%         else
%             %if not bone increment
%             cur_val = cur_val+1;
%         end
%     end
%     
%     %d3-
%     cur_val = bone_dist_maps(6,cs1,cs2,cs3);
%     while(cur_val<r)
%         %test for bone
%         if(padded_bone_map(pcs1,pcs2,pcs3-cur_val))
%             %if bone set val and break
%             %bone_dist_maps(2,cs1,cs2,cs3) = cur_val;
%             %num_bone_directions(cs1,cs2,cs3) = num_bone_directions(cs1,cs2,cs3)+1;
%             %update points along line
%             bone_dist_maps(6,cs1,cs2,cs3-(cur_val-1):cs3) = r+1;
%             nearby_bone_map_by_direction(6,cs1,cs2,cs3-(cur_val-1):cs3) = 1;
%             break;
%         else
%             %if not bone increment
%             cur_val = cur_val+1;
%         end
%     end
    
    
end


for i_d = 1:num_dir
    num_dir_bone_found_in = num_dir_bone_found_in+ squeeze(nearby_bone_map_by_direction(i_d,:,:,:));
end
%imshow3D(num_bone_directions,[0 max(max(max(num_bone_directions)))]);
%imshow3D(squeeze(bone_dir_breakdowns(1,:,:,:)),[0 1]);
internal_map = zeros(size(num_dir_bone_found_in));
internal_map(num_dir_bone_found_in==num_dir) = 1;
updated_bone_map = bone_map;
updated_bone_map(internal_map==1) = 1;


end

function linspace_subs = CustomLinspace(start_point,iterator,end_point, map_size_3d)
%returns indices into the 3D map 
    cur_point = start_point;
    linspace_subs = [];
    linspace_subs = [linspace_subs; cur_point(1),cur_point(2),cur_point(3)];
    while(~isequal(cur_point,end_point))
        cur_point = cur_point+iterator;
        linspace_subs = [linspace_subs; cur_point(1),cur_point(2),cur_point(3)];
    
    end
    
end