function score = ScoreEvaluate2(test_ct3D, dog_filt_1,dog_filt_2,pos_points,neg_points, dog_filt_3, dog_filt_4)
    
    score_adjust = 0;
    if(dog_filt_1<=0)
        dog_filt_1 = 1e-4;
        score_adjust = score_adjust+500;
    end
    if(dog_filt_2<=0)
        dog_filt_2 = 1e-4;
        score_adjust = score_adjust+500;
    end
    if(dog_filt_3<=0)
        dog_filt_3 = 1e-4;
        score_adjust = score_adjust+500;
    end
    if(dog_filt_4<=0)
        dog_filt_4 = 1e-4;
        score_adjust = score_adjust+500;
    end



    %generate dog for region
    gf1 = imgaussfilt3(test_ct3D,dog_filt_1);
    gf2 = imgaussfilt3(test_ct3D,dog_filt_2);
    dog_filt = gf1-gf2;
    
    gf3 = imgaussfilt3(test_ct3D,dog_filt_3);
    gf4 = imgaussfilt3(test_ct3D,dog_filt_4);
    dog_filt2 = gf3-gf4;
    
    
    good_points_scores = 0;
    %pos points
    for i_pp = 1:size(pos_points,1)
        cur_point_ind = sub2ind(size(dog_filt),pos_points(i_pp,1),pos_points(i_pp,2),pos_points(i_pp,3));
        good_points_scores = good_points_scores + (2000+dog_filt(cur_point_ind)) + (2000+dog_filt2(cur_point_ind));
    end
    
    bad_points_scores = 0;
    %pos points
    for i_pp = 1:size(neg_points,1)
        cur_point_ind = sub2ind(size(dog_filt),neg_points(i_pp,1),neg_points(i_pp,2),neg_points(i_pp,3));
       
        bad_points_scores = bad_points_scores + (2000+dog_filt(cur_point_ind)) + (2000+dog_filt2(cur_point_ind));
    end
    
    
    total_score = good_points_scores - bad_points_scores;
    
    score = -1*total_score; %for minimisation
    score = score + score_adjust;
end