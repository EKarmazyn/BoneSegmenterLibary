%first test DB 1
clear all
close all
dbstop if error


for dataset=3:2:10


start_time = tic;
dt = string(Configuration.TorsoDatasets);
dtt = dt(dataset);
SegmentBoneV4(char(dtt), 'pre_volume_classifier_save', 1, 'fixed_initial_threshold',1175, 'dataset', dataset);
tElapsed = toc(start_time);
fprintf("Time-Taken-Overall: " + num2str(tElapsed)+"\n");


end