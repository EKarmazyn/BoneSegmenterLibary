
spec_map = zero_map;

if(use_indices)
    %for i_tv = 1:length(target_vol)
    %    spec_map(target_vol_inds{i_tv}) = 1;
    %end
     spec_map(target_vol_inds) = 1;
else
    
    spec_map(final_labelled_map==target_vol) = 1;
    
end


multi_target = false;
if(length(target_vol)>1)
    multi_target = true;
end

%change bounds to focus on target vol

target_indices = find(spec_map);
[target_sub1, target_sub2, target_sub3] = ind2sub(size(final_labelled_map),target_indices);

border_radius = 40;

min_x1 = max(1, min(target_sub1)-border_radius);
min_x2 = max(1, min(target_sub2)-border_radius);
min_x3 = max(1, min(target_sub3)-border_radius);

max_x1 = min(size(ct3D,1), max(target_sub1)+border_radius);
max_x2 = min(size(ct3D,2), max(target_sub2)+border_radius);
max_x3 = min(size(ct3D,3), max(target_sub3)+border_radius);


spec_map = spec_map(min_x1:max_x1, min_x2:max_x2, min_x3:max_x3);

pos = get(f, 'Position'); %// gives x left, y bottom, width, height
close(f);

if(multi_target)
    f = figure('Name', "multi_target: ");

else
    f = figure('Name', "cur_target_vol: " + num2str(target_vol));

end
set(f, 'Position', pos);
figure(f)

%clf(f)
hold  on
[hpat1] = ViewSeperateVolumes(binary_labelled_map(min_x1:max_x1, min_x2:max_x2, min_x3:max_x3), X_ct, min_x1:max_x1, min_x2:max_x2, min_x3:max_x3);
hpat1.FaceAlpha  = 0;
hpat1.EdgeAlpha  = 0.1;




[hpat2] = ViewSeperateVolumes(spec_map, X_ct, min_x1:max_x1, min_x2:max_x2, min_x3:max_x3);
if(iscell(hpat2))
    hpat2 = hpat2{1};
end
hpat2.FaceColor = 'g';
hpat2.FaceAlpha  = 0.2;
hpat2.EdgeColor = 'g';
hpat2.EdgeAlpha  = 0.2;


% if(multi_target == false)
%     save("C:\Users\mazna\Documents\nl\U\P\Report\FinalReport\FigureData\HardSeperate\IndVolumeAnalysis\" + num2str(target_vol), 'f');
% end
