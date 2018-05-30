%load
clear all
close all
db stop
% 

 dataset = 4;
 dt = Configuration.TorsoDatasets;
 ddt = dt(dataset);
 [connected_volumes_array, X_ct] = SegmentBoneV7(char(ddt),'SaveFolder','C:\Users\mazna\Documents\nl\U\P\Report\FinalReport\FigureData\GOOD\','save_all',true);


%[connected_volumes_arraySK, X_ct2] = SegmentBoneV7(Configuration.ExampleHeadDataPath,'SaveFolder','C:\Users\mazna\Documents\nl\U\P\Report\FinalReport\FigureData\GOOD\','save_all',true);
