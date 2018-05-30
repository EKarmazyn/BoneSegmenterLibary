%load
clear all
close all

dataset = 1;
dt = Configuration.TorsoDatasets;
ddt = dt(dataset);
SegmentBoneV6(char(ddt),'fixed_initial_threshold',1175,  'dataset', 1, 'ss_save', true, 'save_append','_1_test_fT', 'save_result', true);

dataset = 15;
dt = Configuration.TorsoDatasets;
%ddt = dt(dataset);
SegmentBoneV6(Configuration.ExampleHeadDataPath,fixed_initial_threshold',1175,  'dataset', 1, 'ss_save', true, 'save_append','_head_test_fT', 'save_result', true);

