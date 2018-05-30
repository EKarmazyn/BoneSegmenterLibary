function [predictions, scores] = PC3_CustomModel(input_data_table)
%PC3_CUSTOMMODEL custom bone or not bone predictions 
%we expect 1-10 bone connected volumes

%as number of voxels gets large, max intensity should get largeish
for row = 1:size(input_data_table,1)
    cur_row = input_data_table(row,:);
    (5000/cur_row.NumberOfVoxels)
end




end

