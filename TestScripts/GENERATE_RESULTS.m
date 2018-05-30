clear all 
close all
for dataset = 1:10
    %load data
    dt = string(Configuration.TorsoDatasets);
    dtt = dt(dataset);

    SegmentBoneV5(char(dtt), 'dataset', dataset, 'save_results', 1);

end

