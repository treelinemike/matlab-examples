
% dataVectors:
%   m = # rows = # datapoints
%   n = # cols = # of features
%   NO LABELS
function queryLabels = kNearestNeighbors(trainData,trainLabels,queryData,K)

% initialize predicted query label vector
queryLabels = zeros(size(queryData,1),1);

for queryIdx = 1:size(queryData,1)
    
    vectorDiffs = trainData' - repmat(queryData(queryIdx,:)',1,size(trainData,1))
    
    
    
    queryLabels(queryIdx,1) = 0;
end


end


