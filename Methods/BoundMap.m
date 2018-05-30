function [bound_map,ind1,ind2,ind3, full_size] = BoundMap(map,border)
%BOUNDMAP bounds in the imput map (such that zeros are removed, keeps
%border

full_size = size(map);

% d=1;
% i=0;
% while (i<full_size(d))
%     i=i+1;
%     slice =map(i,:,:);
%     %tic
%     if(any(slice(:)))
%         ind1 = i;
%         break;
%     end
%     %toc
% end
% 
% i=full_size(d);
% while (i>0)
%     i=i-1;
%     %tic
%     slice =map(i,:,:);
%     %tic
%     if(any(slice(:)))
%         ind1_b = i;
%         break;
%     end
%     %toc
% end
% 
% %%%%%%%%%%%%%%%%%%55
% d=2;
% i=0;
% while (i<full_size(d))
%     i=i+1;
%     slice =map(:,i,:);
%     %tic
%     if(any(slice(:)))
%         ind2 = i;
%         break;
%     end
% end
% 
% i=full_size(d);
% while (i>0)
%     i=i-1;
%     slice =map(:,i,:);
%     %tic
%     if(any(slice(:)))
%         ind2_b = i;
%         break;
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%
% 
% d=3;
% i=0;
% while (i<full_size(d))
%     i=i+1;
%     slice =map(:,:,i);
%     %tic
%     if(any(slice(:)))
%         ind3 = i;
%         break;
%     end
% end
% 
% i=full_size(d);
% while (i>0)
%     i=i-1;
%      slice =map(:,:,i);
%     %tic
%     if(any(slice(:)))
%         ind3_b = i;
%         break;
%     end
% end

%%% UPDATE TO USE THIS AS FASTER WHEN TIME
% %minimal crop
%     [xx,yy,zz] = ind2sub(size(ct3D),connected_volumes(current_cc).IndVoxelList);
%     x_bot = max([1, min(xx)-1]);
%     x_top = min([size(ct3D,1),max(xx) + 1]);
%     y_bot = max([1, min(yy)-1]);
%     y_top = min([size(ct3D,2),max(yy) + 1]);
%     z_bot = max([1, min(zz)-1]);
%     z_top = min([size(ct3D,3),max(zz) + 1]);



aaa = regionprops(map,'BoundingBox');
if(isempty(aaa))
    bound_map = [];
    ind1 = [];
    ind2 = [];
    ind3 = [];
    full_size = [];
    return
end
ind1 = ceil(aaa.BoundingBox(2));
ind2 = ceil(aaa.BoundingBox(1));
ind3 = ceil(aaa.BoundingBox(3));
ind1_b = ind1+aaa.BoundingBox(5)-1;
ind2_b = ind2+aaa.BoundingBox(4)-1;
ind3_b = ind3+aaa.BoundingBox(6)-1;


ind1 = max(1,ind1-border);
ind2 = max(1,ind2-border);
ind3 = max(1,ind3-border);

ind1_b = min(full_size(1),ind1_b+border);
ind2_b = min(full_size(2),ind2_b+border);
ind3_b = min(full_size(3),ind3_b+border);

bound_map = map(ind1:ind1_b, ind2:ind2_b, ind3:ind3_b);

end

