% function to compute AIC and BIC from a k-Means clustering classification
function [AIC, BIC] = kMeansAICBIC(dataVectors,clusterIDs,meanLocations)

% determine k, n, and a list of unique cluster IDs (this will be trivial)
k = size(meanLocations,1);
n = size(dataVectors,1);
uniqueClusters = unique(clusterIDs);

% initialize RSS to zero
RSS = 0;

% sum RSS from each cluster
for clusterIdx = 1:length(uniqueClusters)    % length(uniqueClusters) should be k
    
    % extract datapoints in cluster
    thisCluster = uniqueClusters(clusterIdx); % uniqueClusters(clusterIdx) == clusterIdx)
    dataInCluster = dataVectors((clusterIDs == thisCluster),:)';
    
    % compute sume of squared differences for all points in cluster
    % relative to cluster mean
    vectorDiff = dataInCluster - repmat(meanLocations(clusterIdx,:)',1,size(dataInCluster,2));
    sumSqDiff = sum(diag(vectorDiff'*vectorDiff));
    
    % add sum of square differences for this cluster to overall RSS
    RSS = RSS  + sumSqDiff;
end
RSS
% compute AIC and BIC from k, n, and RSS
AIC = n*log( RSS / n ) + 2*k;
BIC = n*log( RSS / n ) + k*log( n );
end