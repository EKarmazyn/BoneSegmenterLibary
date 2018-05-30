
dataset = 1;
%load data
dt = string(Configuration.TorsoDatasets);
dtt = dt(dataset);

[~, ct3D, CT_dimension_spacing, X_ct] = LoadDICOM(char(dtt), false);


load("C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\PreSeperation\" + num2str(dataset) + ".mat");

for i_r = 1:3
    act_r = 2*i_r;
    [final_labelled_map{i_r}] = PreClassifySeperate2(ct3D,binary_map_filled, CT_dimension_spacing,X_ct, true, act_r);
end


%load("C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\PostSeperation\" + num2str(dataset) + "_no_coating.mat")

for i = 1:3
    figure
    ViewSeperateVolumes(final_labelled_map{1},X_ct);
end
