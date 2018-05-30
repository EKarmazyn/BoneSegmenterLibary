%bone_points

points{1} = bone_points;


%non_bone_points

points{2} = non_bone_points;

slice = 1;
%z = Slices{1}.Position(3);

x_ct = (CT_dimension_spacing(1)/2):CT_dimension_spacing(1):CT_dimension_spacing(1)*(size(ct3D,1)-(CT_dimension_spacing(1)/2));
y_ct = (CT_dimension_spacing(2)/2):CT_dimension_spacing(2):CT_dimension_spacing(2)*(size(ct3D,2)-(CT_dimension_spacing(2)/2));
z_ct = (CT_dimension_spacing(3)/2):CT_dimension_spacing(3):CT_dimension_spacing(3)*(size(ct3D,3));
X_ct = {x_ct,y_ct,z_ct};

cur_set_positions = [];

cur_points = points{1};
for cp = 1:size(cur_points,1)
    bone_set_positions(cp,:) = [x_ct(cur_points(cp,2)),y_ct(cur_points(cp,1)),z_ct(cur_points(cp,3))];
end

cur_points = points{2};
for cp = 1:size(cur_points,1)
    non_bone_set_positions(cp,:) = [x_ct(cur_points(cp,2)),y_ct(cur_points(cp,1)),z_ct(cur_points(cp,3))];
end
    
[edgeInformationBone] = EdgeInformation(ct3D, X_ct, bone_set_positions ,1,0);
[edgeInformationNonBone] = EdgeInformation(ct3D, X_ct, non_bone_set_positions ,1,0);

bone_dt = EdgeInfo2DataTable(edgeInformationBone,1,1);
non_bone_dt = EdgeInfo2DataTable(edgeInformationNonBone,1,0);

full_table = [bone_dt; non_bone_dt];


