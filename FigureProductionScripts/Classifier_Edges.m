good_edges = load('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\EdgeData\Good\edge_information_10000_spread.mat');

bad_edges = load('C:\Users\mazna\Documents\nl\U\P\Code\MatLab\SavedData\EdgeData\Incorrect\edge_information_10000_spread.mat');



num_total_good_edges = length(good_edges.edge_infomation);
num_total_bad_edges = length(bad_edges.edge_infomation);

to_test =2000;

training_good_edges = good_edges.edge_infomation(1:num_total_good_edges-to_test);
training_bad_edges = bad_edges.edge_infomation(1:num_total_bad_edges-to_test);

test_good_edges = good_edges.edge_infomation(num_total_good_edges-to_test+1:num_total_good_edges);
test_bad_edges = bad_edges.edge_infomation(num_total_bad_edges-to_test+1:num_total_bad_edges);


eval1 = [];
eval2 = [];
eval3 = [];

modelFit = [];

largestGradient = [];

modelSigma = [];

height = [];

width = [];

edgeType = [];


cur_set = training_good_edges;
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
cur_set = training_bad_edges;

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


%str = string(edgeType);

%craft table

input_data_table = table(eval1,eval2,eval3,modelFit,largestGradient,modelSigma,height,width,edgeType);


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
% trainingData = input_data_table;
% inputTable = trainingData;
% predictorNames = {'eval1', 'eval2', 'eval3', 'modelFit', 'largestGradient', 'modelSigma', 'height', 'width'};
% predictors = inputTable(:, predictorNames);
% response = inputTable.edgeType;
% isCategoricalPredictor = [false, false, false, false, false, false, false, false];
% 
% % Train a classifier
% % This code specifies all the classifier options and trains the classifier.
% classificationKNN = fitcknn(...
%     predictors, ...
%     response, ...
%     'Distance', 'Euclidean', ...
%     'Exponent', [], ...
%     'NumNeighbors', 1, ...
%     'DistanceWeight', 'Equal', ...
%     'Standardize', true, ...
%     'ClassNames', [0; 1]);
% 
% % Create the result struct with predict function
% predictorExtractionFcn = @(t) t(:, predictorNames);
% knnPredictFcn = @(x) predict(classificationKNN, x);
% trainedClassifier.predictFcn = @(x) knnPredictFcn(predictorExtractionFcn(x));
% 
% % Add additional fields to the result struct
% trainedClassifier.RequiredVariables = {'eval1', 'eval2', 'eval3', 'modelFit', 'largestGradient', 'modelSigma', 'height', 'width'};
% trainedClassifier.ClassificationKNN = classificationKNN;
% trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2017b.';
% trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');
% 
% % Extract predictors and response
% % This code processes the data into the right shape for training the
% % model.
% inputTable = trainingData;
% predictorNames = {'eval1', 'eval2', 'eval3', 'modelFit', 'largestGradient', 'modelSigma', 'height', 'width'};
% predictors = inputTable(:, predictorNames);
% response = inputTable.edgeType;
% isCategoricalPredictor = [false, false, false, false, false, false, false, false];
% 
% % Perform cross-validation
% partitionedModel = crossval(trainedClassifier.ClassificationKNN, 'KFold', 5);
% 
% % Compute validation predictions
% [validationPredictions, validationScores] = kfoldPredict(partitionedModel);
% 
% % Compute validation accuracy
% validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');




