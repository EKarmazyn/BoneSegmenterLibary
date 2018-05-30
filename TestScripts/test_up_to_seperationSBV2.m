
for dataset = 1:10
    tic
    dt = string(Configuration.TorsoDatasets);
    dtt = dt(dataset);
    [new_labelled_map, X_ct, minimum_bone_intensity] = SegmentBoneV2(char(dtt), 'end_after_seperation',1,'end_after_fuzzy_find_bone',0);
    toc
    save("G:\\COverflowDump\\nl_\\U\\P\\Data\\SeperationOutputs\\\SetsForClassification\\labelled_"+num2str(dataset)+"_"+datestr(now, 'HH-MM-dd-mmm-yyyy'),'new_labelled_map','X_ct');
    
end

%ViewSeperateVolumes(new_labelled_map,X_ct);

close all
dataset = 3
    dt = string(Configuration.TorsoDatasets);
        dtt = dt(dataset);
    [~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(char(dtt), false);
    [valid_mask] = CreateValidAreaMask(ct3D);
    [bone_model_map, minimum_bone_intensity] = FuzzyFindBone(ct3D, valid_mask);
