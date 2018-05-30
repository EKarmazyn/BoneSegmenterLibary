function [edgeInformation] = EdgeInformation(ct3D, X_ct, set_positions ,sample_spacing,display)
%EDGEINFORMATION Extracts information about edges
%%% Inputs:
%       ct3D : Intensity volume
%       X_ct : A cell array containing the tranlation from indices to
%       axis units. (i.e. 1->-150mm, 2-> -148.5mm etc)
%       normals of the edge set
%       positions of the edge set (in actual positions)
%       spacing of samples to be taken from the set (default 1)
%       whether to display information and plot graphs about the edges


%%Deal with optional inputs
if(nargin<4)
    sample_spacing = 1;
    display = false;
elseif(nargin<5)
    display = false;
end



%parameters
lr = 10;
loc_spacing = 1;
r = 8; %11x11x11 with mm spacing
spacing = 0.2;
al_spacing = 0.1;

edge_infomation = [];

%num_samples = 10000;
total_dp = size(set_positions,1);


for dp = 1:sample_spacing:(total_dp - mod(total_dp,sample_spacing))
    %interp to 2r+1 by 2r+1 by 2r+1 cube around the vertex
    
    
    center = set_positions(dp,:);
    x = center(1)-lr*loc_spacing:loc_spacing:center(1)+lr*loc_spacing;
    y = center(2)-lr*loc_spacing:loc_spacing:center(2)+lr*loc_spacing;
    z = center(3)-lr*loc_spacing:loc_spacing:center(3)+lr*loc_spacing;
    
    X = {x,y,z};
    %X_ct = {x_ct,y_ct,z_ct};
    first_x_ct = nan(3,1);
    last_x_ct = nan(3,1);
    for i = 1:3
        %first_x_ct
        shift = X_ct{i}-X{i}(1);
        shift(shift<0) = 0;
        temp = find(shift);
        first_x_ct(i) = temp(1) - 1;

        %last_x_ct
        shift = X_ct{i}-X{i}(end);
        shift(shift<0) = 0;
        temp = find(shift);
        last_x_ct(i) = temp(1);
        
    end
    
    %create x,y,z vector in ct3D co-ords 
    X_local_ct = cell(1,3);
    for i = 1:3
        X_local_ct{i} = first_x_ct(i):last_x_ct(i);
    end
    
    %create xR,yR,zR real locations 
    xR = X_ct{1}(X_local_ct{1});
    yR = X_ct{2}(X_local_ct{2});
    zR = X_ct{3}(X_local_ct{3});
    
    %ct_z_range = ct3D(:,:,X_local_ct{3}(1):X_local_ct{3}(end));
    ct_surrounds = ct3D(X_local_ct{2}(1):X_local_ct{2}(end),X_local_ct{1}(1):X_local_ct{1}(end),X_local_ct{3}(1):X_local_ct{3}(end));
    
    local_gridInterp = griddedInterpolant({xR,yR,zR},ct_surrounds,'makima');
    

    x = center(1)-r*spacing:spacing:center(1)+r*spacing;
    y = center(2)-r*spacing:spacing:center(2)+r*spacing;
    z = center(3)-r*spacing:spacing:center(3)+r*spacing;
    
    %get locality
    [mgx,mgy,mgz] = ndgrid(x,y,z);
    
    locality = local_gridInterp(mgx,mgy,mgz);
    
    %find locality shape
    
    grad_vol_r = 3; %in voxels not length units
    smoothing_sigma = grad_vol_r/2;
    %smoothed_locality = imgaussfilt3(locality,smoothing_sigma);
    
    
     %flip dxdy test
    %center point gradient calculations
    [dy,dx,dz] = gradient(locality);
    
    
    dx2 = dx.*dx;
    dy2 = dy.*dy;
    dz2 = dz.*dz;
    
    %p = mvncdf([mgx(:) mgy(:) mgz(:)], center);
    %p3 = reshape(p,2*r+1,2*r+1,2*r+1);
    %dx = imgaussfilt3(dx,smoothing_sigma, 'Filtersize', 2*r+1,
    weights = NDgaussian(3,r,smoothing_sigma);
    
    
    
    %calculate structure tensor
    %Ix = line_gradient(1);
    %Iy = line_gradient(2);
    %Iz = line_gradient(3);
    sixy = sum(sum(sum(weights.*(dx.*dy))));
    sixz = sum(sum(sum(weights.*(dx.*dz))));
    siyz = sum(sum(sum(weights.*(dy.*dz))));
    
    structure_tensor  = [  sum(sum(sum(weights.*dx2))), sixy, sixz ;
                           sixy, sum(sum(sum(weights.*dy2))) siyz ;
                           sixz, siyz sum(sum(sum(weights.*dz2))) ];
    
    [V,D] = eig(structure_tensor);
    eigVals{dp} = diag(D);
    
    max_grad = V(:,3);
    norm_line_gradient = max_grad;
    
    point_line_gradient = [dx(r+1,r+1,r+1);dy(r+1,r+1,r+1);dz(r+1,r+1,r+1)];
    norm_point_line_gradient = point_line_gradient./(norm(point_line_gradient));
    
    angle = acos(dot(norm_line_gradient,norm_point_line_gradient));
    deg_an = rad2deg(angle);
    %rotate to align with edge
    %unit_normal
    %rotation = vrrotvec(norm_line_gradient,[0,0,1]);
    %rm = vrrotvec2mat(rotation);

    
    
   
    
    %ofset center so we arent on the edge but just outside it!
    
    center = center' - 20*al_spacing*norm_point_line_gradient;
    
    
    pos_step = norm_line_gradient*al_spacing;
    neg_step = norm_line_gradient*(-1*al_spacing);
    num_steps = (lr/2)/al_spacing;
    line_positions = linspaceNDim(center+(num_steps*neg_step), center+(num_steps*pos_step), num_steps*2+1); %center'+(num_steps*neg_step):pos_step:center'+(num_steps*pos_step);
    
  
    current_line = local_gridInterp(line_positions(1,:),line_positions(2,:),line_positions(3,:));
    
    
    pos_step = norm_point_line_gradient*al_spacing;
    neg_step = norm_point_line_gradient*(-1*al_spacing);
    num_steps = (lr/2)/al_spacing;
    line_positions = linspaceNDim(center+(num_steps*neg_step), center+(num_steps*pos_step), num_steps*2+1); %center'+(num_steps*neg_step):pos_step:center'+(num_steps*pos_step);
    
  
    sec_line = local_gridInterp(line_positions(1,:),line_positions(2,:),line_positions(3,:));
    
    
    
    %edge_line_vals(dp,:) = current_line;
    
    
    if(display)
        figure
        plot(current_line,'b');
        hold on;
        plot(sec_line,'r');
        title(num2str(deg_an));
        waitforbuttonpress
        close
        
    end
    current_line = sec_line;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %calc line properties
    %total drop
    %highest_point
    [pks, locs, w, p] = findpeaks(current_line);
    
    %first peak to right (with reasonable p)
    pk = ones(size(pks));
    pk(p<5) = 0;
    pk(locs<num_steps+1) = 0;
    max_val_loc_temp = find(pk);
    if(length(max_val_loc_temp)>0)
        
        max_val_index = max_val_loc_temp(1);
        max_val = pks(max_val_index);
        max_val_loc = locs(max_val_index);
    else
        max_val_loc = length(current_line);
        max_val = current_line(max_val_loc);
    end
    
    
    
    %[max_val, max_val_loc] = max(current_line);
    
    %lowest point (to left of highest)
    [pks, locs, w, p] = findpeaks(current_line*-1);
    pk = ones(size(pks));
    pk(p<1) = 0;
    pk(locs>num_steps+1) = 0;
    min_val_loc_temp = find(pk);
    if(length(min_val_loc_temp)>0)
        
        min_val_index = min_val_loc_temp(end);
        min_val = -1*pks(min_val_index);
        min_val_loc = locs(min_val_index);
    else
        min_val_loc = 1;
        min_val = current_line(min_val_loc);
    end
    %[min_val, min_val_loc] = min(current_line(1:max_val_loc));
    
    %remove bottom edge
    bot_val = (max_val-min_val)*0.05 + min_val;
    temp = current_line;
    temp(1:min_val_loc) = nan;
    [~, bot_val_loc] = min(abs(temp-bot_val));
    bot_val = current_line(bot_val_loc);
    
    if(bot_val_loc < min_val_loc)
        bot_val_loc = min_val_loc;
        bot_val = min_val;
    end
    
    %max grad
    dl = gradient(current_line);
    [max_grad, max_grad_loc] = max(dl(min_val_loc:max_val_loc));
    
    
    %edge sharpness?
    
    
    
    line_segment_double = current_line(min_val_loc:max_val_loc) - min_val;
    %line_segment_double = [line_segment_double flip(line_segment_double)];
    %how like a half normal it is (and then properties)
    %phat = mle(flip(current_line(min_val_loc:max_val_loc)),'distribution','hn');
    %gaussEqn = 'a*exp(-((x-b)/c)^2)';
    %startPoints = [1.5 max_val_loc-bot_val_loc 10];
    %[f1, gof, output] = fit((1:length(line_segment_double))',line_segment_double', gaussEqn, 'Start' , startPoints);
    
    scale = max_val - min_val;
    
    cdf = strcat(num2str(scale/2),'*(1+erf((x-a)/(b*sqrt(2))))');
    mean_guess = (max_val_loc-bot_val_loc)/2+(bot_val_loc - min_val_loc);
    startPoints = [mean_guess 1];
    [f1, gof, output] = fit((1:length(line_segment_double))',line_segment_double', cdf, 'Start' , startPoints);
    
    if(display)
        figure
        plot(f1,(1:length(line_segment_double))',line_segment_double')
        waitforbuttonpress
        close
    end
    %
    
    
    
    %calc variance from hn (euclidean)
    %only counting values between bot and max (not min and max)
    x_ls = (bot_val_loc-min_val_loc)+1:(max_val_loc-min_val_loc)-1;
    residuals = output.residuals(x_ls);
    std_hn = std(residuals);
    
    %structure tensor eVal
    
    %structure tensor eVal ratio
    
    
    
    %store vals
    edge_infomation(dp).EdgeValues = current_line;
    edge_infomation(dp).Max.Value = max_val;
    edge_infomation(dp).Max.Location = max_val_loc;
    edge_infomation(dp).Min.Value = min_val;
    edge_infomation(dp).Min.Location = min_val_loc;
    edge_infomation(dp).Bot.Value = bot_val;
    edge_infomation(dp).Bot.Location = bot_val_loc;
    edge_infomation(dp).LargestGradient.Value = max_grad;
    edge_infomation(dp).LargestGradient.Location = max_grad_loc;
    edge_infomation(dp).Model = f1;
    edge_infomation(dp).ModelResidualSTD = std_hn;
    edge_infomation(dp).EigenValues = diag(D);
    edge_infomation(dp).Position = set_positions(dp,:);
    
    
    
   
    %
    
end


edge_information_long_form = edge_infomation;
edge_infomation = edge_information_long_form(arrayfun(@(s) ~isempty(s.Max),edge_information_long_form));



edgeInformation = edge_infomation;
end

