function [bone_filepath,not_bone_filepath, unknown_filepath] = SaveBone3D(sizeI, connected_volumes_array, axis_vectors)
%SAVEBONE3D saves three ply files, one for bone, one for not bone and one
%for unknown


%first generate the binary maps
bone_map = zeros(sizeI);
not_bone_map = zeros(sizeI);
unknown_map = zeros(sizeI);

for co = 1:length(connected_volumes_array)
    switch(connected_volumes_array(co).Bone)
        case -1
            not_bone_map(connected_volumes_array(co).IndVoxelList) = 1; 
        case 0
            unknown_map(connected_volumes_array(co).IndVoxelList) = 1; 
        case 1
            bone_map(connected_volumes_array(co).IndVoxelList) = 1; 
        otherwise
            error('Bone properpty not set in a connected volume object!')
    end
end


iso = isosurface(axis_vectors{1}, axis_vectors{2}, axis_vectors{3},not_bone_map,0.9);
not_bone_filepath = strcat('C:\Users\mazna\Documents\nl\U\P\Data\GeneratedPly\FINAL\','not_bone_map',datestr(now, 'HH-MM-dd-mmm-yyyy'),'.ply');
plywrite(not_bone_filepath,iso.faces,iso.vertices);

iso = isosurface(axis_vectors{1}, axis_vectors{2}, axis_vectors{3},bone_map,0.9);
bone_filepath = strcat('C:\Users\mazna\Documents\nl\U\P\Data\GeneratedPly\FINAL\','bone_map',datestr(now, 'HH-MM-dd-mmm-yyyy'),'.ply');
plywrite(bone_filepath,iso.faces,iso.vertices);

iso = isosurface(axis_vectors{1}, axis_vectors{2}, axis_vectors{3},unknown_map,0.9);
unknown_filepath = strcat('C:\Users\mazna\Documents\nl\U\P\Data\GeneratedPly\FINAL\','unknown_map',datestr(now, 'HH-MM-dd-mmm-yyyy'),'.ply');
plywrite(unknown_filepath,iso.faces,iso.vertices);

end

