

addpath('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\PotentiallyUsefulLibraries\NIfTI_20140122');
addpath('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\PotentiallyUsefulLibraries');

%%Load DICOM (intensity data)
[Slices, ct3D, CT_dimension_spacing] = LoadDICOM(Configuration.ExampleDicomDataPath, false);


x_ct = (CT_dimension_spacing(1)/2):CT_dimension_spacing(1):CT_dimension_spacing(1)*(size(ct3D,1)-(CT_dimension_spacing(1)/2));
y_ct = (CT_dimension_spacing(2)/2):CT_dimension_spacing(2):CT_dimension_spacing(2)*(size(ct3D,2)-(CT_dimension_spacing(2)/2));
z_ct = (CT_dimension_spacing(3)/2):CT_dimension_spacing(3):CT_dimension_spacing(3)*(size(ct3D,3));
%X = {x,y,z};
X_ct = {x_ct,y_ct,z_ct};

%Load .nii bone model
bone_set = false;
display = false;

if(bone_set)
    nii = load_nii(Configuration.NifitSkeletonDataPath);
else
    %Load .nii incorrect heart-ish model
    nii = load_nii('C:\Users\mazna\Documents\nl\U\P\Data\Example\dicom_data\example_incorrect_volume\example_heart.nii');
end

%



bone_point_cloud_dense = nii.img;
bone_point_cloud_dense(bone_point_cloud_dense ~= 0) = 1;
%FLIP to align orientation to ct3D
bone_point_cloud_dense = flip(bone_point_cloud_dense,3);
bone_point_cloud_dense = imrotate(bone_point_cloud_dense, 90);
bone_point_cloud_dense = flip(bone_point_cloud_dense,2);


%isosurface
iso = isosurface(bone_point_cloud_dense, 0.9);
ct_coord_vertices = iso.vertices;
%figure
%patch(iso);

vertices = nan(size(ct_coord_vertices));

%transform vertices to correct co-ordinates
for dp = 1:size(ct_coord_vertices,1)
    for i = 1:3
        vertices(dp,i) = X_ct{i}(round(ct_coord_vertices(dp,i)));
    end
    
   
end

%first calculate normal directions
n = isonormals(x_ct, y_ct, z_ct, ct3D, vertices);

%now magnitudes
nm = cellfun(@norm, num2cell(n, 2));


if(bone_set)
    %remove small normals <30
   [sni] = find(nm > 50);   %THIS IS NORMAL OPERATION
    
    %%%%% TESTING / LOOKING AT bad edges
     %[sni] = find(nm < 50);

    good_set_magnitudes = nm(sni);
    good_set_positions = vertices(sni,:);
    good_set_normals = n(sni,:);
    
    
else
    good_set_magnitudes = nm;
    good_set_positions = vertices;
    good_set_normals = n;
end

%remove vertices mear the edge of the ct data;
lr = 10; %this is the radius we need to leave clear


bounds = nan(3,2);
for i = 1:3
    bounds(i,1) = lr;
    bounds(i,2) = X_ct{i}(end)-lr;
    
end


edge_point_is = [];
for dp = 1:size(good_set_positions,1)
    if(    good_set_positions(dp,1) < bounds(1,1) || good_set_positions(dp,1) > bounds(1,2)...
        || good_set_positions(dp,2) < bounds(2,1) || good_set_positions(dp,2) > bounds(2,2)...
        || good_set_positions(dp,3) < bounds(3,1) || good_set_positions(dp,3) > bounds(3,2))
        edge_point_is = [edge_point_is dp];
    end
end


%remove edge points
good_set_magnitudes(edge_point_is) = [];
good_set_positions(edge_point_is,:) = [];
good_set_normals(edge_point_is,:) = [];


%x_ct = (CT_dimension_spacing(1)/2):CT_dimension_spacing(1):CT_dimension_spacing(1)*(size(ct3D,1)-(CT_dimension_spacing(1)/2));
%y_ct = (CT_dimension_spacing(2)/2):CT_dimension_spacing(2):CT_dimension_spacing(2)*(size(ct3D,2)-(CT_dimension_spacing(2)/2));
%z_ct = (CT_dimension_spacing(3)/2):CT_dimension_spacing(3):CT_dimension_spacing(3)*(size(ct3D,3));
gridVectCells = {x_ct,y_ct,z_ct};

gridInterp = griddedInterpolant(gridVectCells, ct3D, 'makima');

%dataset to test:
%good set
%note that the positions given in good_set_positions are given in a unit
%co-ordinate assumption->


num_samples = 10000;
total_dp = size(good_set_normals,1);
sample_spacing = (total_dp - mod(total_dp,num_samples))/num_samples;


bad_edge_info = EdgeInformation(ct3D,X_ct,good_set_positions,sample_spacing);

%good_edges = good_edge_info;
bad_edges = bad_edge_info;
%bad_edges = load('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\EdgeData\Incorrect\edge_information_10000_spread.mat');

num_total_good_edges = length(good_edges);
num_total_bad_edges = length(bad_edges);

to_test =2000;

training_good_edges = good_edges(1:num_total_good_edges-to_test);
training_bad_edges = bad_edges(1:num_total_bad_edges-to_test);

test_good_edges = good_edges(num_total_good_edges-to_test+1:num_total_good_edges);
test_bad_edges = bad_edges(num_total_bad_edges-to_test+1:num_total_bad_edges);

tge_dt = EdgeInfo2DataTable(training_good_edges,1,1);
tbe_dt = EdgeInfo2DataTable(training_bad_edges,1,0);

input_data_table = [tge_dt;tbe_dt];


eval1 = [];
eval2 = [];
eval3 = [];

modelFit = [];

largestGradient = [];

modelSigma = [];

height = [];

width = [];

edgeType = [];

cur_set = test_good_edges;


l = length(cur_set);


for i = 1:l
    eval1(i,1) = cur_set(i).EigenValues(1);
    eval2(i,1) = cur_set(i).EigenValues(2);
    eval3(i,1) = cur_set(i).EigenValues(3);
    
    modelFit(i,1) = cur_set(i).ModelResidualSTD;
    
    largestGradient(i,1) = cur_set(i).LargestGradient.Value;
    
    modelSigma(i,1) = cur_set(i).Model.b;
    
    height(i,1) = cur_set(i).Max.Value - cur_set(i).Min.Value;
    
    width(i,1) = cur_set(i).Max.Location - cur_set(i).Min.Location;
    
    edgeType(i,1) = 1;
end

cur_set = [];
cur_set = test_bad_edges;


l2 = length(cur_set);


for j = 1:l2
    i = j+l;
    eval1(i,1) = cur_set(j).EigenValues(1);
    eval2(i,1) = cur_set(j).EigenValues(2);
    eval3(i,1) = cur_set(j).EigenValues(3);
    
    modelFit(i,1) = cur_set(j).ModelResidualSTD;
    
    largestGradient(i,1) = cur_set(j).LargestGradient.Value;
    
    modelSigma(i,1) = cur_set(j).Model.b;
    
    height(i,1) = cur_set(j).Max.Value - cur_set(j).Min.Value;
    
    width(i,1) = cur_set(j).Max.Location - cur_set(j).Min.Location;
    
    edgeType(i,1) = 0;
end


test_data_table = table(eval1,eval2,eval3,modelFit,largestGradient,modelSigma,height,width,edgeType);


%%Train model here and export trainedmodel


yfit=trainedModel.predictFcn(test_data_table);
should_be_one = yfit(1:length(test_good_edges));
correct_ratio_for_bone = sum(should_be_one) / length(should_be_one);

should_be_zero = yfit(length(test_good_edges)+1:end);
correct_ratio_for_heart = (length(test_bad_edges)-sum(should_be_zero)) / length(should_be_zero);

