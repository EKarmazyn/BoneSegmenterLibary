%load .mat

dataset = 2;
dt = string(Configuration.TorsoDatasets);
dtt = dt(dataset);

[~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(char(dtt), false);


load("C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\PreSeperation\" + num2str(dataset) + ".mat");

save_path = "C:\Users\mazna\Documents\nl\U\P\Data\GeneratedPly\pre_classify_vols\";

%seperate step
xy_thresh = 15
z_thresh = 0
[connected_volumes_array, fully_labelled] = PreClassificationSeperate(ct3D,binary_map_filled, X_ct, CT_dimension_spacing, xy_thresh,z_thresh);


%save all and then some, with and without smoothing

min_large = 500;

zero_map = uint8(zeros(size(fully_labelled)));
binary_bone_map = zero_map;
binary_bone_map(fully_labelled>0) = 1;



%full not smoothed
tic
iso = isosurface(X_ct{1}, X_ct{2}, X_ct{3},binary_bone_map,0.95);
fp = save_path + "fullNotSmoothed95_"+ datestr(now, 'HH-MM-dd-mmm-yyyy') + ".ply";
plywrite(fp,iso.faces,iso.vertices);
toc

tic
iso = isosurface(X_ct{1}, X_ct{2}, X_ct{3},binary_bone_map,0.5);
fp = save_path + "fullNotSmoothed50_"+ datestr(now, 'HH-MM-dd-mmm-yyyy') + ".ply";
plywrite(fp,iso.faces,iso.vertices);
toc

%only large & seperated not smoothed method 1(simple)
to_save_list = [];
for i_cv = 1:length(connected_volumes_array)
    if(connected_volumes_array(i_cv).NumVoxels > min_large)
        to_save_list = [to_save_list; i_cv];
    end
end


fprintf("saving all large cvs seperately"+"\n");
TT = tic;


x1= X_ct{1};
x2 = X_ct{2};
x3 = X_ct{3};

parfor i_to_save = 1:length(to_save_list)
    
    fprintf("individual cv save: "+ num2str(i_to_save)+"\n");
    T = tic;
    
    %fprintf("generating zero_map_copy");
    %T = tic;
    cur_map = zero_map;
    %t = toc(T);
    %fprintf(".        Time-Taken: " + num2str(t)+"\n");
    
    %fprintf("filling zero_map_copy");
    %T = tic;
    cur_map(connected_volumes_array(i_to_save).IndVoxelList) = 1;
    %t = toc(T);
    %fprintf(".        Time-Taken: " + num2str(t)+"\n");
    
   % fprintf("iso surf 0.95");
    %T = tic;
    iso = isosurface(x1, x2, x3,cur_map,0.95);
    %t = toc(T);
    %fprintf(".        Time-Taken: " + num2str(t)+"\n");
    
    
    
    fp = save_path + "large500_" + num2str(i_to_save) +"_largeSeperatedNotSmoothed95_"  + datestr(now, 'HH-MM-dd-mmm-yyyy') + ".ply";
    
    
    
    
    %fprintf("plywrite");
    %T = tic;
    plywrite(fp,iso.faces,iso.vertices);
    %t = toc(T);
    %fprintf(".        Time-Taken: " + num2str(t)+"\n");
    
    
    t = toc(T);
    fprintf(".        Time-Taken: " + num2str(t)+"\n");
    
end
t = toc(TT);
fprintf(".        Time-Taken: " + num2str(t)+"\n");
    
    
    
%only large & seperated not smoothed method 2(more complex below)
%%% first smooth and threshold, 
%%% then CC analysis to check for large vols
%%% next identify cvs which are in these vols and add to list
%%% save these cvs

%full smoothed



%only large & seperated smoothed



% iso = isosurface(axis_vectors{1}, axis_vectors{2}, axis_vectors{3},not_bone_map,0.999);
% not_bone_filepath = strcat('C:\Users\mazna\Documents\nl\U\P\Data\GeneratedPly\','not_bone_map',datestr(now, 'HH-MM-dd-mmm-yyyy'),'.ply');
% plywrite(not_bone_filepath,iso.faces,iso.vertices);
