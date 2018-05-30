function [internal_map] = RemoteInternalFill(input_mask, bone_map,search_r,allowed_misses)
%REMOTEINTERNALFILL Summary of this function goes here
% 

asm = NET.addAssembly('C:\Users\mazna\Documents\nl\U\P\Code\CS\dllForMatlab\dllForMatlab\bin\x64\Release\dllForMatlab.dll');

tic
int_fill = dllForMatlab.InternalFill(NET.convertArray(search_r,'System.Byte'),NET.convertArray(allowed_misses,'System.Byte'),NET.convertArray(size(bone_map),'System.Int32'));
toc
r = search_r;
padded_bone_map = ones(size(bone_map,1)+2*r,size(bone_map,2)+2*r,size(bone_map,3)+2*r);
padded_bone_map(r+1:size(bone_map,1)+r,r+1:size(bone_map,2)+r,r+1:size(bone_map,3)+r) = bone_map;
padded_bone_map=uint8(padded_bone_map);

in_arr = find(input_mask==1);
[s1, s2, s3] = ind2sub(size(bone_map), in_arr);



%v = NET.convertArray(in_arr,'System.Byte',[3*3*3]);
m1 = NET.convertArray(padded_bone_map,'System.Byte');
%m2 = NET.convertArray(in_arr,'System.Int32');
s1 = NET.convertArray(s1,'System.Int32');
s2 = NET.convertArray(s2,'System.Int32');
s3 = NET.convertArray(s3,'System.Int32');
%dirs_found = int_fill.Run(m1,s1,s2,s3);
tic
ids_net = int_fill.RunToIds(m1,s1,s2,s3);
toc
ids_list = int32(ids_net);
internal_map= uint8(zeros(size(input_mask)));
for i = 1:length(ids_list)%+1 for ml indexing!
    internal_map(in_arr(ids_list(i)+1)) = 1;
end










%internal_map = uint8(dirs_found)
% a=1;
%  tic
%  temp = GetRange(dirs_found_arr,0,dirs_found_arr.Count);
%  dArr = ToArray(temp);
%  num_dir_bone_found_in= uint8(zeros(size(input_mask)));
%  cell_ar = {};
%  for i_d = 1:26
%     cell_ar{i_d} = uint8(dArr(i_d));
% end
% for i_d = 1:26
%     num_dir_bone_found_in = num_dir_bone_found_in+ squeeze(cell_ar{i_d}(:,:,:));
% end
% count_map = uint8(zeros(size(input_mask)));
% for i=1:dirs_found_arr.Count
%     count_map = count_map + uint8(dArr(i));
% end
% toc
% %SaveForTransfer
% delete('C:\Users\mazna\Documents\nl\U\P\Code\CS\TransferFolder\finished_writing');
% fileID = fopen('C:\Users\mazna\Documents\nl\U\P\Code\CS\TransferFolder\config.txt','w');
% fprintf(fileID,'%s',num2str(search_r) + " " + num2str(allowed_misses) +" \n");
% fclose(fileID);
% 
% fileID = fopen('C:\Users\mazna\Documents\nl\U\P\Code\CS\TransferFolder\input_mask','w');
% 
% fwrite(fileID,input_mask)
% fclose(fileID);
% 
% fileID = fopen('C:\Users\mazna\Documents\nl\U\P\Code\CS\TransferFolder\bone_map','w');
% 
% %fwrite(fileID,bone_map)
% %fclose(fileID);
% 
% save('C:\Users\mazna\Documents\nl\U\P\Code\CS\TransferFolder\input_mask.mat','input_mask');
% save('C:\Users\mazna\Documents\nl\U\P\Code\CS\TransferFolder\bone_map.mat','bone_map');
% pause(0.25);
% 
% %call exe
% 
% 
% 
% 
% while(~exist('plot.m', 'file') )
%     pause(1);
% end
% 
% 
% %GrabFromTransfer
% load('C:\Users\mazna\Documents\nl\U\P\Code\CS\TransferFolder\upd_bone_map.mat');
% 
% 
% %load upd_bone_map
% 
% 
% internal_map = upd_bone_map-bone_map;
end

