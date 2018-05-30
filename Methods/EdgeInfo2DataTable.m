function [data_table] = EdgeInfo2DataTable(edge_information, edge_type, abs_edge_type)
%EDGEINFO2DATATABLE Summary of this function goes here
cur_set = [];
cur_set = edge_information;

unknown = false;
if(nargin < 2)
    unknown = true;

end

if(nargin>2)
    single_type = true;
end


l2 = length(cur_set);


for j = 1:l2
    %i = j+l;
    eval1(j,1) = cur_set(j).EigenValues(1);
    eval2(j,1) = cur_set(j).EigenValues(2);
    eval3(j,1) = cur_set(j).EigenValues(3);
    
    modelFit(j,1) = cur_set(j).ModelResidualSTD;
    
    largestGradient(j,1) = cur_set(j).LargestGradient.Value;
    
    modelSigma(j,1) = cur_set(j).Model.b;
    
    height(j,1) = cur_set(j).Max.Value - cur_set(j).Min.Value;
    
    width(j,1) = cur_set(j).Max.Location - cur_set(j).Min.Location;
    
    top(j,1) = cur_set(j).Max.Value;
    
    pos1(j,1) = cur_set(j).Position(1);
    pos2(j,1) = cur_set(j).Position(2);
    pos3(j,1) = cur_set(j).Position(3);
    
    if(~unknown)
        if(~single_type)
            edgeType(j,1) = cur_set(j).EdgeType;
        else
            edgeType(j,1) = abs_edge_type;
        end
    end
    
end

if(unknown)
    
    data_table = table(eval1,eval2,eval3,modelFit,largestGradient,modelSigma,height,width, top, pos1,pos2,pos3);    
else
    
    data_table = table(eval1,eval2,eval3,modelFit,largestGradient,modelSigma,height,width, top, edgeType, pos1,pos2,pos3);
end


end

