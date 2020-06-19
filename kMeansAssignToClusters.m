function clusterIDs = kMeansAssignToClusters(dataVectors,meanLocations)

clusterIDs = zeros( size(dataVectors,1), 1);

for dataIdx = 1:size(dataVectors,1)
    thisDataLoc = dataVectors(dataIdx,:)';
    vectorDiff = repmat(thisDataLoc,1,size(meanLocations,1)) - meanLocations';
    [~,minIdx] = min(sqrt(diag(vectorDiff'*vectorDiff)));  % don't really need sqrt() here
    clusterIDs(dataIdx,1) = minIdx;
end

end
