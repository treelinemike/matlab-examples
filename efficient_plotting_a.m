% There are lots of ways to improve plotting efficieny in MATLAB to display
% data in close to real time (assuming reasonably slow update rates). This
% demo shows one way of updating data directly in the plot handle which can
% be quicker than re- or over-plotting pre-existing data. Note that the X
% and Y axes are NOT autoscaled and we pre-allocate storage for all data we
% expect to recieve. Alternatively this could be implemented with a rolling
% window similar to an oscilloscope display.

% restart
close all; clear all; clc;

% plot after recieving how many datapoints?
cyclesBetweenPlotting = 5;

% generate simulated data
t = 0:0.01:5;
y = sin(2*pi*t);

% initialize plot
figure;
hold on; grid on;
xlim([0,5]);
ylim([-1.5 1.5]);

% plot null data using pre-allocated NaN vectors
ph = plot(nan(1,1000),nan(1,1000),'-','LineWidth',1.6,'Color',[0.0 0.0 0.8]);

% simulate recieving data 
dataIdx = 1;
readBufferT = zeros(1,1000);
readBufferY = zeros(1,1000);
readBufferIdx = 0;
plotOffset = 0;

% step through simulated input data stream
while(dataIdx <= length(t))
    
    % read one datapoint (byte? word?) from the input stream
    readBufferIdx = readBufferIdx + 1;
    readBufferT(readBufferIdx) = t(dataIdx);
    readBufferY(readBufferIdx) = y(dataIdx);
    
    % update plotting on a less frequent basis
    if(~mod(dataIdx,cyclesBetweenPlotting))
        
        % update data in plot handle and redraw
        ph.XData(plotOffset+(1:readBufferIdx)) = readBufferT(1:readBufferIdx);
        ph.YData(plotOffset+(1:readBufferIdx)) = readBufferY(1:readBufferIdx);
        drawnow;
        
        % update index variables
        plotOffset = plotOffset + readBufferIdx;
        readBufferIdx = 0;
    end
    
    % increment to next datapoint
    dataIdx = dataIdx + 1;
    pause(0.01);  % simulate sampling at ~100Hz (not accurate! to do this correctly we'd wait here for an already-running timer to expire)
end
