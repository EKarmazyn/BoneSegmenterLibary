function  ViewBone3D(labelled_map, connected_volumes_array,X_ct)
%VIEWBONE3D labelled map as from bwlabeln,
% connected volumes in array need to have the Bone value set (to -1,0 or 1)

%clone
view_map = labelled_map;
view_map(labelled_map == 0) = nan;

%generate color map
cmap = jet(3);

for co = 1:length(connected_volumes_array)
    switch(connected_volumes_array(co).Bone)
        case -1
            view_map(connected_volumes_array(co).IndVoxelList) = 1; %red
        case 0
            view_map(connected_volumes_array(co).IndVoxelList) = 3; %blue
        case 1
            view_map(connected_volumes_array(co).IndVoxelList) = 2; %green
        otherwise
            error('Bone properpty not set in a connected volume object!')
            
    end
    
    
end
ViewSeperateVolumes(view_map,X_ct);
%PATCH_3Darray(view_map);


end

