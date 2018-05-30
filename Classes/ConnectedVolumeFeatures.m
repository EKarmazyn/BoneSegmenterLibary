classdef ConnectedVolumeFeatures < Feature
    %CONNECTEDVOLUMEFEATURE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        MaxIntensity
        MeanIntensity
        MedianIntensity
        IntensityVariance
        Intensity025
        Intensity25
        Intensity75
        Intensity975
        DistanceOfCoMToDefBoneSeedCO
        DistanceOfCoMToProbBoneSeedCO
        VoxelPCA_Data
        EdgeInfomation %array
        PercentageEdgeClassifierBH
        WeightedPercentageEdgeRegressionBH
        PercentageEdgeClassifierNew
        BasicPercentageSurfaceEdge
        MinDistToOtherCO
        MinDistToOtherCO_Diameter_Adjusted
        NumberOfVoxels
        ErodeFraction
        
        EquivDiameter
        Extent
        EVects
        EVals
        ConvexVolume
        Solidity
        SurfaceArea
        
        
        BoneSeed %0,1,2 for none, prob, def
        
        
    end
    
    methods
        
    end
end

