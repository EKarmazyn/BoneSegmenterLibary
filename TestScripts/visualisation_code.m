
image = seperated_map;
image(basic_surfaces==1) = 2;
figure
imshow3D(image,[0 max(max(max(image)))]);

image = cur_surf_points;
figure
imshow3D(image,[0 max(max(max(image)))]);

figure
ViewSeperateVolumes(labelled_new_model, X_ct)


image = ct3D;
figure
imshow3D(image,[0 max(max(max(image)))]);