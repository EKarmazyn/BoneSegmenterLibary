function [Slice, data3D, dimension_spacing, X_ct] = LoadDICOM(dicom_data_directoryPath, primitive_correction)
%LOADDICOM Summary of this function goes here
%   Detailed explanation goes here

% dimension_spacing contains the unit spacing 
% data3D

%% Setup

%clear all;
rtnFolder = pwd;

%% Load Slices

directoryPath_exampleDicomData = dicom_data_directoryPath;
cd(directoryPath_exampleDicomData);
% Directory contains many slices
slicesFI = dir('*.dcm');
if(length(slicesFI) == 0)
    slicesFI = dir('IMG*');
end
% slices appear to be sorted by name already, may need to add functionality
% to ensure this.

num_slices = length(slicesFI);

Slice = [];


slice_arbitary_origin = 0;
slice_arbitary_axis_pointer = slice_arbitary_origin;

for i_slice = 1:num_slices
    Slice(i_slice).Info = dicominfo(slicesFI(i_slice).name);
    Slice(i_slice).Data = dicomread(Slice(i_slice).Info);
    %%Slice(i_slice).Data(Slice(i_slice).Data == 0) = NaN; % maybe find
    %%circle and make results outside circle NaN.
    % position of the top left corner of the data (image). Also, 
    % ImageOrientationPatient contains the orientation.
    Slice(i_slice).Position = Slice(i_slice).Info.ImagePositionPatient;
end


% Get general info (for all the slices, assuming they have the same info)
HighBit = Slice(1).Info.HighBit;
LargestPosiblePixelValue = 2^(HighBit+1)-1;
PixelSpacing = Slice(1).Info.PixelSpacing;
Rows = Slice(1).Info.Rows;
Columns = Slice(1).Info.Columns;
NominalSliceThickness = Slice(1).Info.SliceThickness;
% i think this 1.5 is 1.5mm, the slices therefore overlap, which could be 
% important to note.
SpacingBetweenSlices = abs(Slice(2).Position(3) - Slice(1).Position(3)); 
TotalLength = abs(Slice(end).Position(end) - Slice(1).Position(end));
% a method to go from overlapping slices to re-constructed non-overlapping
% slices might be useful, half the data from one slice is also in the slice
% adjacent on that slice, and vice versa. 



cd(rtnFolder);

%% Primitive correction
if(nargin>1)
    if(primitive_correction)
        %line up histograms!
        
%         
%         
%         for i_slice = 1:num_slices
%             xraysum(i_slice) = sum(sum(Slice(i_slice).Data));
%         end
%         
%         slice_avg_sum = sum(xraysum)/num_slices;
%         
%         
%         
%         for i_slice = 1:num_slices
%             adj_factor = slice_avg_sum / xraysum(i_slice);
%             for x = 1:Columns
%                 for y = 1:Rows
%                     Slice(i_slice).Data(x,y) = round(Slice(i_slice).Data(x,y)*adj_factor);
%                 end
%             end
%         end
    end
end

%% Store data as 3d matrix

data3D = nan(Rows, Columns, num_slices);
for i_slice = 1:num_slices
    data3D(:,:,i_slice) = Slice(i_slice).Data;
end
data3D = uint16(data3D);

%%
dimension_spacing = [PixelSpacing(1), PixelSpacing(2), SpacingBetweenSlices];

x_ct = (dimension_spacing(1)/2):dimension_spacing(1):dimension_spacing(1)*(size(data3D,1)-(dimension_spacing(1)/2));
y_ct = (dimension_spacing(2)/2):dimension_spacing(2):dimension_spacing(2)*(size(data3D,2)-(dimension_spacing(2)/2));
z_ct = (dimension_spacing(3)/2):dimension_spacing(3):dimension_spacing(3)*(size(data3D,3));
%X = {x,y,z};
X_ct = {x_ct,y_ct,z_ct};

end

