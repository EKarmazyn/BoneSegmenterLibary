function   UserMarkBone(fuzzy_bone_map,ct3D, connected_volumes)
%USERMARKBONE 



%plot sparse surface points to check
%figure
%hold on;
%imshow3D(valid_cc);
% Create push button
 
color_cc_map = zeros(size(fuzzy_bone_map,1),size(fuzzy_bone_map,2),size(fuzzy_bone_map,3),3);

for color = 1:3
    
    switch(color)
        case 1
            bone_tag = -1;
        case 2
            bone_tag = 1;
        case 3
            bone_tag = 0;
    end
    
    all_points_ind_list = ConnectedVolume.GetAllPoints(connected_volumes, bone_tag);
    
    %ind 2 sub
    [co_sub{1},co_sub{2},co_sub{3}] = ind2sub(size(fuzzy_bone_map),all_points_ind_list);


    esub_add{4} = squeeze(zeros(length(co_sub{1}),1) + color); %blue, +1 is red, +2 green
    if(size(esub_add{4},1) == 0)
        esub_add{4}=[];
    end
    %sub to ind

    co_ind_add = sub2ind(size(color_cc_map),co_sub{1},co_sub{2},co_sub{3},esub_add{4});

    %update color map

    color_cc_map(co_ind_add) = 1;
end







figure;
ref_axes_handle = imshow3D(ct3D);
%user mark volumes as bone
%cur_fig = figure('units','normalized','outerposition',[0 0 1 1])
cur_fig  = figure;
main_axes_handle = imshow3D(color_cc_map, [], true, @call_back_fcn, ref_axes_handle, ct3D);
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);


%toggle button for bone or not bone

tb = uicontrol(gcf, 'Style', 'togglebutton', 'String', 'Bone', 'Position', [30 400 200 100]);

%secondary window for viewing ct data


figure(cur_fig);
%wrap slider callback function to also update slice of other figure




while(true)
    try
        [Y,X,~,~,S] = ginputc(100,'Color','r');
    

        %update ;bone; 
        i_ind = [];
        for p = 1:length(Y)
            i_ind(p) = sub2ind(size(fuzzy_bone_map),round(X(p)),round(Y(p)),S(p));
        end
        %i_ind = 
    catch
        continue
    end
    
    for p = 1:length(i_ind)
        i_co = ConnectedVolume.GetConnectedVolumeInd(connected_volumes, i_ind(p));

        if(get(tb,'Value'))
             bone_val = 1;
             color_val = 2;
        else
            bone_val = -1;
            color_val = 1;
        end

        if(i_co ~= -1)
            connected_volumes(i_co).Bone = bone_val;

            %ind 2 sub
            [co_sub{1},co_sub{2},co_sub{3}] = ind2sub(size(ct3D),connected_volumes(i_co).IndVoxelList);

            %expand sub
            % for D = 1:3
            %     esub{D} = [co_sub{D}; co_sub{D}; co_sub{D}];
            % end
            esub_remove{4} = zeros(length(co_sub{1}),1) +3; %remove unknown blue
            esub_add{4} = zeros(length(co_sub{1}),1) + color_val; %add color representing (green bone, red not bone)

            %sub to ind
            co_ind_rem = sub2ind(size(color_cc_map),co_sub{1},co_sub{2},co_sub{3},esub_remove{4});
            co_ind_add = sub2ind(size(color_cc_map),co_sub{1},co_sub{2},co_sub{3},esub_add{4});

            %update color map
            color_cc_map(co_ind_rem) = 0;
            color_cc_map(co_ind_add) = 1;
        end
    end


    cla(main_axes_handle) ;

    main_axis_handle = imshow3D(color_cc_map, [], true, @call_back_fcn, ref_axes_handle, ct3D);
    %set(get(ref_axes_handle,'children'),'cdata',squeeze(ct3D(:,:,S,:)));
    %set(get(main_axis_handle,'children'),'cdata',squeeze(color_cc_map(:,:,S,:)));
    figure(cur_fig);
    
end

end


function call_back_fcn(S, ref_axis, refImg)

        set(get(ref_axis,'children'),'cdata',squeeze(refImg(:,:,S,:)));
        set(gcf,'UserData',S);
%         caxis([Rmin Rmax])
%         if sno > 1
%             set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
%         else
%             set(stxthand, 'String', '2D image');
%         end
end
