classdef ConnectedVolume < handle
    %CONNECTEDVOLUME 
    
    properties
        SurfPoints
        IndVoxelList
        IndSurfList
        NumVoxels
        BinaryMap
        Bone %-1 = not bone, 1 = bone, 0 = not known
        BoneSeed %2 = def bone seed, 1 = prob bone seed, 0 = unknown, -1 = no
        Features %obj handle
        Centroid
        EdgeInfo %on the SurfPoints
        VolumeSize %0 unknown, 1-5, v_small to v_large -1 too small
    end
    
    methods
        
        function obj = ConnectedVolume(IndVoxelList)
            %CONNECTEDVOLUME Construct an instance of this class
            %   Detailed explanation goes here
            if(nargin>0)
                obj.IndVoxelList = IndVoxelList;
                obj.NumVoxels = length(IndVoxelList);
                obj.Bone = 0;%unknown
                obj.BoneSeed = 0;
            end
        end
        
        function GenerateEdgeInfo(obj, ct3D, X_ct)
            obj.EdgeInfo = EdgeInformation(ct3D, X_ct, obj.SurfPoints);
        end
        
        function GenerateIndividualFeatures(obj, ct3D, BH_class_model, BH_reg_model)
            if(isempty(obj.Features))
                cvf = ConnectedVolumeFeatures;
                %Features = cvf;
                obj.Features = cvf;
            else
                cvf = obj.Features;
            end
            
            
            
            intensities = ct3D(obj.IndVoxelList);
            y = quantile(double(intensities), [0.025 0.25 0.5 0.75 0.975]);
            % MaxIntensity
            cvf.MaxIntensity = max(intensities);
            % MeanIntensity
            cvf.MeanIntensity = mean(intensities);
            % MedianIntensity
            cvf.MedianIntensity = median(intensities);
            % IntensityVariance
            cvf.IntensityVariance = var(double(intensities));
            % Intensity025
            cvf.Intensity025 = y(1);
            % Intensity25
            cvf.Intensity25 = y(2);
            % Intensity75
            cvf.Intensity75 = y(4);
            % Intensity975
            cvf.Intensity975 = y(5);
            
            % NumberOfVoxels
            cvf.NumberOfVoxels = obj.NumVoxels;
            
            % BoneSeed
            cvf.BoneSeed = obj.BoneSeed;
            
            %GenerateRegionPropsFeatures VoxelPCA_Data
            
            
            % GROUP FEATURE DistanceOfCoMToDefBoneSeedCO
            % GROUP FEATURE DistanceOfCoMToProbBoneSeedCO
            
            
            % EdgeInfomation %array NOT DONE
            % PercentageEdgeClassifierBH NOT DONE
            % WeightedPercentageEdgeRegressionBH NOT DONE
            % GROUP FEATURE PercentageEdgeClassifierNew
            % GROUP FEATURE MinDistToOtherCO
            
        end
        
        function [x,y,z] = CenterOfMass(obj)
            error('not implemented');
        end
        
        
        
        
        
        function [binary_map] = ToBinaryMap(obj, szVol)
            if(~isnan(obj.BinaryMap))
                binary_map = obj.BinaryMap;
                return
            end
            binary_map = zeros(szVol);
            [x,y,z] = ind2sub(szVol,obj.IndVoxelList);
            binary_map(x,y,z) = 1;
            obj.BinaryMap = binary_map;
        end
        
        
      
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        function GenerateBasicPercentageSurfEdge(obj)
            
        end
    end
    
    methods(Static)
        
        function [labelled_map] = ToLabelledMap(labelled_map)
            connected_volumes_array(1,max(labelled_map(:))) = ConnectedVolume;
            for i = 1:max(labelled_map(:))
                connected_volumes_array(1,i) = ConnectedVolume(find(labelled_map==i));
            end
            
        end
        
        function [connected_volumes_array] = FromLabelledMap(labelled_map)
            connected_volumes_array(1,max(labelled_map(:))) = ConnectedVolume;
            for i = 1:max(labelled_map(:))
                connected_volumes_array(1,i) = ConnectedVolume(find(labelled_map==i));
            end
            
        end
        
        function [connected_volumes_array] = FromLabelledIndices(cell_array_of_ind_lists)
            connected_volumes_array(1,length(cell_array_of_ind_lists)) = ConnectedVolume;
            for i = 1:length(cell_array_of_ind_lists)
                connected_volumes_array(1,i) = ConnectedVolume(cell_array_of_ind_lists{i});
            end
            
        end
        
        function GenerateAllEdgeInfo(connected_volumes_array, ct3D, X_ct)
            for co = 1:length(connected_volumes_array)
                cur_co = connected_volumes_array(co);
                cur_co.GenerateEdgeInfo(ct3D,X_ct);
            end
        end
        
        function [full_table, bone_table, non_bone_table, labelled_table, unlabelled_table] = ConstructDataTable(connected_volumes_array)
            
            % MaxIntensity
            % MeanIntensity
            % MedianIntensity
            % IntensityVariance
            % Intensity025
            % Intensity25
            % Intensity75
            % Intensity975
            % DistanceOfCoMToDefBoneSeedCO
            % DistanceOfCoMToProbBoneSeedCO
            % VoxelPCA_Data
            % EdgeInfomation %array
            % PercentageEdgeClassifierBH
            % WeightedPercentageEdgeRegressionBH
            % PercentageEdgeClassifierNew
            % MinDistToOtherCO
            % NumberOfVoxels
            % BoneSeed %0,1,2 for none, prob, def
            % BONETYPE %-1,0,1 for not, unknown, bone
            
            for co = 1:length(connected_volumes_array)
                cur_co = connected_volumes_array(co);
                cur_features = cur_co.Features;
                MaxIntensity(co,1) = cur_features.MaxIntensity;
                MeanIntensity(co,1) = cur_features.MeanIntensity;
                MedianIntensity(co,1) = cur_features.MedianIntensity;
                IntensityVariance(co,1) = cur_features.IntensityVariance;
                Intensity025(co,1) = cur_features.Intensity025;
                Intensity25(co,1) = cur_features.Intensity25;
                Intensity75(co,1) = cur_features.Intensity75;
                Intensity975(co,1) = cur_features.Intensity975;
                BoneSeed(co,1) = cur_features.BoneSeed;
                NumberOfVoxels(co,1) = cur_features.NumberOfVoxels;
                ErodeFraction(co,1) = cur_features.ErodeFraction;
                
                EquivDiameter(co,1) = cur_features.EquivDiameter;
                Extent(co,1) = cur_features.Extent;
                %EVect1(co,1) = cur_features.EVects(1);
                %EVect1(co,1) = cur_features.EVects(1);
                %EVect1(co,1) = cur_features.EVects(1);
                EVal1(co,1) = cur_features.EVals{1}(1);
                EVal2(co,1) = cur_features.EVals{1}(2);
                EVal3(co,1) = cur_features.EVals{1}(3);
                ConvexVolume(co,1) = cur_features.ConvexVolume;
                Solidity(co,1) = cur_features.Solidity;
                SurfaceArea(co,1) = cur_features.SurfaceArea;
                PCA1(co,1) = cur_features.VoxelPCA_Data(1);
                PCA2(co,1) = cur_features.VoxelPCA_Data(2);
                PCA3(co,1) = cur_features.VoxelPCA_Data(3);
                
                if(isempty(cur_features.MinDistToOtherCO))
                    cur_features.MinDistToOtherCO = 0;
                end
                
                if(isempty(cur_features.MinDistToOtherCO_Diameter_Adjusted))
                    cur_features.MinDistToOtherCO_Diameter_Adjusted = 0;
                end
                
                cur_features.MinDistToOtherCO_Diameter_Adjusted
                DistToBone(co,1) = cur_features.MinDistToOtherCO;
                DistToBoneADJ(co,1) = cur_features.MinDistToOtherCO_Diameter_Adjusted;
                
                BoneType(co,1) = cur_co.Bone;
                
            end
            
            full_table = table(MaxIntensity,...
                MeanIntensity,...
                MedianIntensity,...
                IntensityVariance,...
                Intensity025,...
                Intensity25,...
                Intensity75,...
                Intensity975,...
                BoneSeed,...
                NumberOfVoxels,...
                EquivDiameter,...
                Extent,...
                EVal1,...
                EVal2,...
                EVal3,...
                Solidity,...
                ConvexVolume,...
                SurfaceArea,...
                PCA1,...
                PCA2,...
                PCA3,...
                DistToBone, ...
                DistToBoneADJ,...
                ErodeFraction,...
                BoneType);
                
            
            
            rows = full_table.BoneType == 1;
            bone_table = full_table(rows, :);
                                
            rows2 = full_table.BoneType == -1;
            non_bone_table = full_table(rows2, :);
            abcd = rows|rows2;
            labelled_table = full_table(rows|rows2, :);
            ind_list = (rows|rows2)+1;
            ind_list = logical(mod(ind_list,2));
            unlabelled_table = full_table(ind_list, :);
        end
        
        function GenerateErodeFraction(connected_volumes_array, binary_map)
            
            eroded = imerode(binary_map, [1 1 1; 1 1 1; 1 1 1]);
            eroded = binary_map - eroded;
            %erodeded_points = find(binary_map-eroded);
            for co = 1:length(connected_volumes_array)
                cur_co = connected_volumes_array(co);
                cur_features = cur_co.Features;
                
                num_eroded = sum(eroded(cur_co.IndVoxelList));
                
                cur_features.ErodeFraction = num_eroded/cur_co.NumVoxels;
                cur_co.Features = cur_features;
                connected_volumes_array(co) = cur_co;
            end
        end
        
        function GenerateAllBasicPercentageSurfEdges(connected_volumes_array)
            for co = 1:length(connected_volumes_array)
                cur_co = connected_volumes_array(co);
                cur_co.GenerateBasicPercentageSurfEdge();
            end
        end
        
        function GenerateRegionPropsFeatures(connected_volumes_array, labelled_binary_map)
            T= tic;
            stats = regionprops3(labelled_binary_map, {"EquivDiameter",...
                                                        "Extent",...
                                                        "PrincipalAxisLength",...
                                                        "EigenVectors",...
                                                        "EigenValues",...
                                                        "ConvexVolume",...
                                                        "Solidity",...
                                                        "SurfaceArea"});
            t = toc(T);
            fprintf("regionprops3  Time-Taken: " + num2str(t)+"\n");



            for co = 1:length(connected_volumes_array)
                T= tic;
                cur_co = connected_volumes_array(co);
                cur_features = cur_co.Features;
                cur_stats = stats(co,:);
                %cur_co.Centroid = cur_stats.Centroid;
                cur_features.EquivDiameter = cur_stats.EquivDiameter;
                cur_features.Extent = cur_stats.Extent;
                cur_features.VoxelPCA_Data = cur_stats.PrincipalAxisLength;
                cur_features.EVects = cur_stats.EigenVectors;
                cur_features.EVals = cur_stats.EigenValues;
                cur_features.ConvexVolume = cur_stats.ConvexVolume;
                cur_features.Solidity = cur_stats.Solidity;
                cur_features.SurfaceArea = cur_stats.SurfaceArea;
                cur_co.Features = cur_features;
                connected_volumes_array(co) = cur_co;
                t = toc(T);
                fprintf("regionprops3-singleCV  Time-Taken: " + num2str(t)+"\n");



            end
        end
        
        function GenerateSpatialFeatures(connected_volumes_array, labelled_map,X_ct )
            
            vol_lengths = zeros(length(connected_volumes_array),1);

            
            for i_label = 1:length(connected_volumes_array)
                vol_lengths(i_label) = connected_volumes_array(i_label).NumVoxels;
            end

            
            %short_vol_min = 237.8906*0.3965; %100 for 1

            %short_label_min = round(short_vol_min/SCALE.vol);
            short_labels = find(vol_lengths<100);

             %rp  = regionprops3(final_labelled_map, {'Centroid','EquivDiameter'});
             %CoM = rp(:,'Centroid');
             %EqDia = rp(:,'EquivDiameter');
            %setup co-ords transform
            %gi = griddedInterpolant(1:size(ct3D,1), 1:size(ct3D,2), 1:size(ct3D,3), ndgrid(X_ct{1},X_ct{2},X_ct{3}));
            gi = {};
            for i = 1:3
                gi{i} = griddedInterpolant(1:size(labelled_map,i), X_ct{i});

            end

            all_indices = ConnectedVolume.GetAllPoints(connected_volumes_array);


            small_indices_arr = {};
            small_ind_labels = [];

            for i_short = 1:length(short_labels)
                cur_label = short_labels(i_short);
                    small_indices_arr{i_short} = connected_volumes_array(cur_label).IndVoxelList;
                    small_ind_labels = [small_ind_labels; ones( connected_volumes_array(cur_label).NumVoxels,1)...
                                            *cur_label];
            end

            %remove cur label indices from all indices list
            all_small_indices = cell2mat(small_indices_arr');
            other_indices = setdiff(all_indices, all_small_indices);

            %use NN tagging to find distance to closest part of the model

            %unlabelled_point_indices = setdiff(unlabelled_point_indices,all_seed_indexes);

            %[ul_arr] = CoM{short_labels, 'Centroid'};
            [ulp_1,ulp_2,ulp_3] = ind2sub(size(labelled_map),all_small_indices);
            %labelled_point_indices = find(labelled_map);
            %labelled_point_indices = all_seed_indexes;
            [lp_1,lp_2,lp_3] = ind2sub(size(labelled_map),other_indices);
            %convert to real positions using x_ct
            ulp_1 = gi{1}(ulp_1);%swapped on purpose!
            ulp_2 = gi{2}(ulp_2);%
            ulp_3 = gi{3}(ulp_3);

            lp_1 = X_ct{1}(lp_1);
            lp_2 = X_ct{2}(lp_2);
            lp_3 = X_ct{3}(lp_3);



            [~, D] = knnsearch([lp_1',lp_2',lp_3'],[ulp_1,ulp_2,ulp_3]);


            non_short_labels = setdiff(1:length(connected_volumes_array),short_labels);
            for i_cv = 1:length(non_short_labels)
                cur_label = non_short_labels(i_cv);
                cur_co = connected_volumes_array(cur_label);
                cur_features = cur_co.Features;
               
                cur_features.MinDistToOtherCO = 0;
                cur_features.MinDistToOtherCO_Diameter_Adjusted = 0;
                cur_co.Features = cur_features;
                connected_volumes_array(cur_label) = cur_co;
            end
            
            for i_cv = 1:length(short_labels)
                
                
                min_a = min(D(small_ind_labels==i_cv));
                if(length(min_a)>1)
                    min_a = min_a(1);
                end
                cur_label = short_labels(i_cv);
                cur_co = connected_volumes_array(cur_label);
                cur_features = cur_co.Features;
               
                cur_features.MinDistToOtherCO = min_a;
                cur_features.MinDistToOtherCO_Diameter_Adjusted = min_a-cur_features.EquivDiameter;
                cur_co.Features = cur_features;
                connected_volumes_array(cur_label) = cur_co;
            end
            
            
        end
        
        function GenerateGroupFeatures(connected_volumes_array)
            
            
            %first get centroids for each connected_volume
            
            centroids = [];
            
            % GROUP FEATURE DistanceOfCoMToDefBoneSeedCO
            % GROUP FEATURE DistanceOfCoMToProbBoneSeedCO
            % GROUP FEATURE PercentageEdgeClassifierNew
            % GROUP FEATURE MinDistToOtherCO
            
            error('notimplemented');
        end
        
        
        function points_ind_list = GetAllPoints(connected_volumes_array, requested_bone_type)
            points_ind_list = [];
            for co = 1:length(connected_volumes_array)
                cur_co = connected_volumes_array(co);
                if(nargin>1)
                    if(cur_co.Bone == requested_bone_type)
                        points_ind_list = [points_ind_list; cur_co.IndVoxelList];
                    end
                else

                    points_ind_list = [points_ind_list; cur_co.IndVoxelList];
                end
                
                
            end
        end
        
      
        
        function AssignSurfacePoints(connected_volumes_array, surface_points_original_indexing, szVol, actual_surface_points)
            %surface points list in index values
            %for each connected volume find all surface points which are
            %members
            ind_surf_points_orig = sub2ind(szVol,surface_points_original_indexing(:,1),surface_points_original_indexing(:,2),surface_points_original_indexing(:,3));
            %tic
            for co = 1:length(connected_volumes_array)
                cur_co = connected_volumes_array(co);
                logical_map = ismember(ind_surf_points_orig, cur_co.IndVoxelList);
                potential_surface_points = actual_surface_points(logical_map, :);
                if(size(potential_surface_points,1) < 1000)
                    cur_co.SurfPoints = potential_surface_points;
                else
                    %take first 1000, they are already randomised
                    cur_co.SurfPoints = potential_surface_points(1:1000,:);
                end
                
                ind_surf_points_orig(logical_map) = [];
            end
            %toc
            
        end
        
        function MarkWithBoneSeeds(connected_volumes_array,prob_bone_seed_binary_map, def_bone_seed_binary_map)
            [~,prob_bone_seed_points] = find(prob_bone_seed_binary_map);
            [~,def_bone_seed_points] = find(def_bone_seed_binary_map);
            for co = 1:length(connected_volumes_array)
                cur_co = connected_volumes_array(co);
                if(any(ismember(def_bone_seed_points, cur_co.IndVoxelList)))
                    cur_co.BoneSeed = 2;
                elseif(any(ismember(prob_bone_seed_points, cur_co.IndVoxelList)))
                    cur_co.BoneSeed = 1;
                else
                    cur_co.BoneSeed = -1;
                end
            end
        end
        
        function i_co = GetConnectedVolumeInd(connected_volumes, i_ind)
            %finds the connected volume  containing the voxel linearly
            %indexed by i_ind
            for it_co = 1:length(connected_volumes)
                if(ismember(i_ind, connected_volumes(it_co).IndVoxelList))
                    i_co = it_co;
                    return;
                end
            end
            i_co = -1;
        end
    end
end

