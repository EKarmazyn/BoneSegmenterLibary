clear all 
close all
 dataset = 1
    %load data
    dt = string(Configuration.TorsoDatasets);
    dtt = dt(dataset);

    SegmentBoneV6(char(dtt), 'dataset', dataset, 'ss_save', true);



