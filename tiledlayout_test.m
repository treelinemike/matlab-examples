% potentially life-changing replacement for subplot?
% added in R2019b

% restart
close all; clear; clc;

% options
numRows = 2;
numCols = 3;

% create a tiled layout
figure;
set(gcf,'Position',[1.898000e+02 1.602000e+02 7.312000e+02 5.248000e+02]);  % adjust as necessary
t = tiledlayout(numRows,numCols);

% add individual tiles
for rowIdx = 1:numRows
    for colIdx = 1:numCols
        thisTileIdx = (rowIdx-1)*numCols+colIdx;
        nexttile(thisTileIdx);
        hold on; grid on;
        spy(eye(thisTileIdx),50);
    end
end

% get rid of that pesky whitespace!
t.Padding = 'none';  % 'normal', 'compact', or 'none'
t.TileSpacing = 'none';  % 'normal', 'compact', or 'none'
