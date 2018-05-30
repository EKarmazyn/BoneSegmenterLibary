%%% Purpose of this script is to produce the classifiers from the training
%%% data in: G:\COverflowDump\nl_\U\P\Data\SeperationOutputs\SetsForClassification
clear all

load('G:\COverflowDump\nl_\U\P\Data\SeperationOutputs\SetsForClassification\labelled_3_02-41-05-Mar-2018.mat')
training_map = new_labelled_map;
X_ct_training = X_ct;

load('G:\COverflowDump\nl_\U\P\Data\SeperationOutputs\SetsForClassification\labelled_4_02-47-05-Mar-2018.mat')
testing_map = new_labelled_map;
X_ct_testing = X_ct;

clear new_labelled_map;

connected_volumes = ConnectedVolume.FromLabelledMap(training_map);

%test old classifier

dataset = 3;
dt = string(Configuration.TorsoDatasets);
dtt = dt(dataset);
[~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(char(dtt), false);%find probable and definate bone seeds
[prob_bone_seed_map, prob_minimum_bone_seed_intensity, def_bone_seed_map, def_bone_seed_min_intensity] = FindBoneSeeds(ct3D, ones(size(ct3D)));

%mark connected volumes with bone seeds
ConnectedVolume.MarkWithBoneSeeds(connected_volumes,prob_bone_seed_map, def_bone_seed_map);


%extract features of connected volumes
for cv = 1:length(connected_volumes)
    connected_volumes(cv).GenerateIndividualFeatures(ct3D);
    
end

ConnectedVolume.GenerateRegionPropsFeatures(connected_volumes, training_map);


%load model
load('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\Models\PC3\new_classi_model2GOOD.mat');


%temp_ind = [1:11, 13:length(connected_volumes)];
[full_table, bone_table, non_bone_table] = ConnectedVolume.ConstructDataTable(connected_volumes);
%input_table = [bone_table; non_bone_table];

%test model
yfit = new2_classi_model.predictFcn(full_table);
tested_connected_volumes = connected_volumes;
%Hightlight3D(fuzzy_bone_map,connected_volumesALT(3),[],[]);
bone_ind = (yfit==1);
not_bone_ind = (yfit==-1);

bone_vols = tested_connected_volumes(bone_ind);
for bv = 1:length(bone_vols)
    bone_vols(bv).Bone = 1;
end

not_bone_vols = tested_connected_volumes(not_bone_ind);
for bv = 1:length(not_bone_vols)
    not_bone_vols(bv).Bone = -1;
end


bone_map = training_map~=0;

Highlight3D(bone_map,tested_connected_volumes(not_bone_ind),tested_connected_volumes(bone_ind),[]);

UserMarkBone(bone_map,ct3D, connected_volumes);

%THE BED IS 20 AND NEEDS REMOVING
connected_volumes(20).Bone=0;%unknown as bed
[full_table, bone_table, non_bone_table] = ConnectedVolume.ConstructDataTable(connected_volumes);


%test model
yfit = trainedModel.predictFcn(full_table);
tested_connected_volumes = connected_volumes;
%Hightlight3D(fuzzy_bone_map,connected_volumesALT(3),[],[]);
bone_ind = (yfit==1);
not_bone_ind = (yfit==-1);

bone_vols = tested_connected_volumes(bone_ind);
for bv = 1:length(bone_vols)
    bone_vols(bv).Bone = 1;
end

not_bone_vols = tested_connected_volumes(not_bone_ind);
for bv = 1:length(not_bone_vols)
    not_bone_vols(bv).Bone = -1;
end


Highlight3D(bone_map,tested_connected_volumes(not_bone_ind),tested_connected_volumes(bone_ind),[]);

%%% TESTING SET

connected_volumesTest = ConnectedVolume.FromLabelledMap(testing_map);


dataset = 4;
dt = string(Configuration.TorsoDatasets);
dtt = dt(dataset);
[~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(char(dtt), false);%find probable and definate bone seeds
[prob_bone_seed_map, prob_minimum_bone_seed_intensity, def_bone_seed_map, def_bone_seed_min_intensity] = FindBoneSeeds(ct3D, ones(size(ct3D)));

%mark connected volumes with bone seeds
ConnectedVolume.MarkWithBoneSeeds(connected_volumesTest,prob_bone_seed_map, def_bone_seed_map);


%extract features of connected volumes
for cv = 1:length(connected_volumesTest)
    connected_volumesTest(cv).GenerateIndividualFeatures(ct3D);
    
end

ConnectedVolume.GenerateRegionPropsFeatures(connected_volumesTest, testing_map);

%temp_ind = [1:11, 13:length(connected_volumes)];
[full_table, bone_table, non_bone_table] = ConnectedVolume.ConstructDataTable(connected_volumesTest);
%input_table = [bone_table; non_bone_table];

%test model
yfit = trainedModel.predictFcn(full_table);
tested_connected_volumes = connected_volumesTest;
%Hightlight3D(fuzzy_bone_map,connected_volumesALT(3),[],[]);
bone_ind = (yfit==1);
not_bone_ind = (yfit==-1);

bone_vols = tested_connected_volumes(bone_ind);
for bv = 1:length(bone_vols)
    bone_vols(bv).Bone = 1;
end

not_bone_vols = tested_connected_volumes(not_bone_ind);
for bv = 1:length(not_bone_vols)
    not_bone_vols(bv).Bone = -1;
end


bone_map = training_map~=0;

Highlight3D(bone_map,tested_connected_volumes(not_bone_ind),tested_connected_volumes(bone_ind),[]);

