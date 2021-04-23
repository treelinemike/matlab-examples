% potentially life-changing replacement for subplot?
% added in R2019b; options revised in 

% restart
close all; clear; clc;

% options
numRows = 2;
numCols = 3;

% create a tiled layout
figure;
set(gcf,'Position',[1.898000e+02 1.602000e+02 7.312000e+02 5.248000e+02]);  % adjust as necessary
t = tiledlayout(numRows,numCols);

% get rid of that pesky whitespace!
t.Padding = 'loose';  % 'loose', 'compact', or 'tight'
t.TileSpacing = 'tight';  % 'loose', 'compact', 'tight', or 'none'


% add individual tiles
ax = [];
for rowIdx = 1:numRows
    for colIdx = 1:numCols
        thisTileIdx = (rowIdx-1)*numCols+colIdx;
        ax(end+1) = nexttile(thisTileIdx);
        hold on; grid on;
        spy(eye(thisTileIdx),50);
        if(rowIdx == 1)
            set(ax(end),'XTickLabel',{});
            xlabel('');
        else
            xlabel('\bfX Axis Label');
        end
        if(colIdx == 1)
            ylabel('\bfY Axis Label');
        else
            set(ax(end),'YTickLabel',{});
        end
    end
end

% link axes
linkaxes(ax,'xy');

% add titles
title(t,'\bfSuper Title');
xlabel(t,'\bfSuper X Axis Label');
ylabel(t,'\bfSuper Y Axis Label');