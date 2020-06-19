
% dataVectors:
%   m = # rows = # datapoints
%   n = # cols = # of features
%   NO LABELS
function [meanLocations, clusterIDs] = kMeansBuildClusters(dataVectors,K,convergenceThreshold)

% flag signaling convergence
runClustering = 1;

% plotting options, etc.
doAnimate = 1;
doSavePlotFrames = 1;
plotFrameIdx = 1;
colors = [.9 0 0; 0 0 .9; 0 .65 0; .75 0 .75; 0 .8 .8; .8 .6 0; 1 .3 0];
colorWhite = [1 1 1];
colorFadeFactor = 0.7;

% find range of data
mSamples  = size(dataVectors,1);
nFeatures = size(dataVectors,2);

% initialize K means to random locations chosen uniformly across range of data
meanLocations = unifrnd( repmat(min(dataVectors),K,1), repmat(max(dataVectors),K,1),K,nFeatures);

% assign each datapoint to closest mean
clusterIDs = kMeansAssignToClusters(dataVectors,meanLocations);

% create a new figure if animation is to be shown
if(doAnimate)
    figure;
end

while( runClustering )
    
    % plot data and cluster mean locations
    if(doAnimate)
        hold off;
        cla;
        hold on; grid on;
        for clusterIdx = 1:K
            thisColorIdx = mod(clusterIdx-1,size(colors,1))+1;
            thisColor = colors(thisColorIdx,:);
            thisColor = thisColor + colorFadeFactor*(colorWhite-thisColor);
            pointsInCluster = dataVectors(find(clusterIDs == clusterIdx),:);
            plot(pointsInCluster(:,1),pointsInCluster(:,2),'.','Color',thisColor,'MarkerSize',15);
        end
        for meanIdx = 1:size(meanLocations,1)
            thisColorIdx = mod(meanIdx-1,size(colors,1))+1;
            thisColor = colors(thisColorIdx,:);
            plot(meanLocations(meanIdx,1),meanLocations(meanIdx,2),'x','MarkerSize',15,'LineWidth',4,'Color',thisColor);
        end
        title(sprintf('2D Visualization of K-Means Clustering w/ K = %d',K));
        if(doSavePlotFrames)
            saveas(gcf,sprintf('frame%03d.png',plotFrameIdx));
            plotFrameIdx = plotFrameIdx + 1;
        else
            pause(0.1);
        end
    end
    
    % reset candidate mean locations and cluster IDs
    newMeanLocations = zeros(size(meanLocations));
    newClusterIDs = zeros(mSamples,1);
    
    % compute candidate mean location updates
    clusterIDKey = unique(clusterIDs);
    for clusterIdx = 1:length(clusterIDKey)
        pointsInCluster = dataVectors((clusterIDs == clusterIDKey(clusterIdx)),:);
        newMeanLocations(clusterIdx,:) = sum(pointsInCluster,1)/size(pointsInCluster,1);
    end
    
    
    % compute Euclidian distance between old and new mean locations
    vectorDiff = meanLocations' - newMeanLocations';
    maxMeanDelta = max(sqrt(diag(vectorDiff'*vectorDiff)));
    
    % does this new clustering produce a significant change in any of the mean locations?
    % if so, update mean locations and data point cluster assignments
    % if not, end algorithm using previous mean locations and data point cluster assignments
    if(maxMeanDelta > convergenceThreshold)
        meanLocations = newMeanLocations;
        clusterIDs = kMeansAssignToClusters(dataVectors,meanLocations);
        
        % if any cluster doesn't have any members, bump it to a random
        % location to try again
        emptySets = setdiff((1:K)',clusterIDs);
        while( ~isempty(emptySets) )
            for emptySetIdx = 1:length(emptySets)
                meanLocations(emptySets(emptySetIdx),:) =  unifrnd( min(dataVectors), max(dataVectors));
            end
            clusterIDs = kMeansAssignToClusters(dataVectors,meanLocations);
            emptySets = setdiff((1:K)',clusterIDs);
        end
        
    else
        runClustering = 0;
    end
    
end  % while( runClustering )

end


