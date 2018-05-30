function [new_connected_volumes_list] = SeperateVolumes(connected_volumes,ct3D)
%SEPERATEVOLUMES voids edges to attempt to seperate volumes by their
%internal surfaces. Returns a more seperated list of connected volumes
r = 2;
SE = strel('sphere',r);

new_connected_volumes_list = [];
for i_cc = 1:length(connected_volumes)
    current_cc = i_cc;
    %extract gray level map
    grey_map = zeros(size(ct3D));
    grey_map(connected_volumes(current_cc).IndVoxelList) = ct3D(connected_volumes(current_cc).IndVoxelList);
    
    
    
    %minimal crop
    [xx,yy,zz] = ind2sub(size(ct3D),connected_volumes(current_cc).IndVoxelList);
    x_bot = max([1, min(xx)-1]);
    x_top = min([size(ct3D,1),max(xx) + 1]);
    y_bot = max([1, min(yy)-1]);
    y_top = min([size(ct3D,2),max(yy) + 1]);
    z_bot = max([1, min(zz)-1]);
    z_top = min([size(ct3D,3),max(zz) + 1]);
    
    grey_map_crop = grey_map(x_bot:x_top, y_bot:y_top, z_bot:z_top);
    
    %extract binary surface
    bin_map = ones(size(grey_map_crop));
    bin_map(grey_map_crop==0) = 0;
    
    surf_map = bin_map - imerode(bin_map, ones(3,3,3));
    connected_volumes(current_cc).IndSurfList = find(surf_map);
    
    %%gradient map
    [gx,gy,gz] = gradient(grey_map_crop);
    gmag2 = gx.^2 + gy.^2 + gz.^2; %%squareroot unnecessary for flat threshold
    
    
    
    %remove near-edge gradients
    %apply mask
    gmag2test = zeros(size(gmag2));
    fbm = find(bin_map);
    gmag2test(fbm) = gmag2(fbm);
    
    %gmag2test(connected_volumes(current_cc).IndSurfList) = 0 ;
    
    


    t1_test = 5e5;

    %threshold
    threshold_grads = zeros(size(gmag2test));


    threshold_grads(gmag2test>t1_test) = 1;

    %view_map = grey_map_crop;
    %view_map(threshold_grads==1) = 3000;
    %figure
    %imshow3D(view_map);

    %%dilate-Erode-dilate
    threshold_grads2 = imdilate(threshold_grads, SE);
    %figure
    %imshow3D(threshold_grads2, [0 max(max(max(threshold_grads)))]);

    %%Check Seperation
    bin_map2 = bin_map;
    bin_map2((threshold_grads2)==1) = 0;

    
    [CC] = bwconncomp(bin_map2);
    %[labelled_binary_mapCC,N] = bwlabeln(bin_map2);
    labelled_binary_mapCC = labelmatrix(CC);
    
    if(isempty(CC.PixelIdxList))
        new_connected_volumes_list = [new_connected_volumes_list; connected_volumes(current_cc)];
        continue;
    end
    
    
    %start_cc
%     start_cc_map = zeros(size(labelled_binary_mapCC));
%     start_cc_map(labelled_binary_mapCC == current_cc) = 1;
%     figure
%     imshow3D(start_cc_map);
    
    
    %unlabelled 
    unlabelled_map = bin_map;
    unlabelled_map(labelled_binary_mapCC~=0) = 0;
    
    %spread till all of original bin_map is full
    full_labelled_binary_map_newCC = labelled_binary_mapCC;
    
    
    %iterative largest first
    %TODO
    cust_dil_radius = 1;
    while(any(unlabelled_map(:)))
        unlabelled_count = sum(unlabelled_map(:))
        for n = 1:length(CC.PixelIdxList)
            spec_map = zeros(size(full_labelled_binary_map_newCC));
            spec_map(full_labelled_binary_map_newCC==n) = 1;
            
            %dil1 = imdilate(spec_map, ones(ln,0,0));
            %dil2 = imdilate(spec_map, ones(0,ln,0));
            %dil3 = imdilate(spec_map, ones(0,0,ln));
            dil = CustomDilate(spec_map, cust_dil_radius);
            upd = dil;% & dil;
            upd(unlabelled_map==0) = 0;
            
            full_labelled_binary_map_newCC(upd==1) = n;
            %update unlabelled_map
            unlabelled_map(upd==1) = 0;

        end
        if(unlabelled_count == sum(unlabelled_map(:)))
            cust_dil_radius = cust_dil_radius + 1;
        end
        
    end
    
    %add new connected volumes to new list
    
    %%stick labelled_binary_map into full_size
    full_size_labelled_map = zeros(size(ct3D));
    full_size_labelled_map(x_bot:x_top, y_bot:y_top, z_bot:z_top) = full_labelled_binary_map_newCC;
    
    for n = 1:length(CC.PixelIdxList)
        new_connected_volumes_list = [new_connected_volumes_list; ConnectedVolume(find(full_size_labelled_map==n))];
    end
end










end

