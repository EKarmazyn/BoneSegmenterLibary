function [] = SaveKeyPointsImgs(kp,sl,label_mat, fp)
%SAVEKEYPOINTSIMGS Summary of this function goes here



side_length = sl;


%kp(2).x = [96];% 186];
%kp(2).y = [217];% 305] ;
%kp(2).z = 291 ;


    
        cur_img = label_mat(kp.x:kp.x+side_length, kp.y:kp.y+side_length, kp.z);
        cur_labels = unique(cur_img);
        con = length(cur_labels)-1;
        labeled = zeros(size(cur_img));
        for i = 1:con
            labeled(cur_img == cur_labels(i+1))= i;
        end
        %f = figure
%         if(max(cur_img(:))>1)
%             cur_img = label2rgb(cur_img,'spring','c','shuffle'); 
%             %imshow(cur_img);
%         else
%             cur_img = logical(cur_img);
%         end
        %imshow(cur_img, [0, max(cur_img(:))]);
        
        %saveas(gcf,"C:\Users\mazna\Documents\nl\U\P\Report\FinalReport\FigureData\kp\im_kp" + num2str(i_kp) + "_s" +  num2str(i) + "_d" +num2str(i) + ".png");
        fn = fp + datestr(now, 'HH-MM-dd-mmm-yyyy') + ".png";
        
        RGB = label2rgb(labeled,parula);
        
        
        imwrite(RGB, char(fn));
        %colormap(f, 'jet');
   
    


end

