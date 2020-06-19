% restart
close all; clearvars; clc;

load kmeansdata
x = X(:,3:4);

meanLocs = kMeansBuildClusters(x,7,eps);
