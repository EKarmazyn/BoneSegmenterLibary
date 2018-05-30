function [] = SaveResults(variable,results_name, figure_handle )
%SAVERESULTS wrapper, Saves results to results folder

save("C:\Users\mazna\Documents\nl\U\P\Data\ML_results\" + results_name, 'variable'); 

if(nargin>2)
    save("C:\Users\mazna\Documents\nl\U\P\Data\ML_results\fig_" + results_name, 'figure_handle'); 

end
end

