% a set of continuous functions can have a discontinuous limit...
% same chart as on cover of real analysis book

% restart
close all; clear all; clc;

nvals = logspace(0,10,30);
xvals = 0:0.001:1;

figure;
hold on; grid on;
xlim([0 1]);
ylim([0 1]);

for nIdx = 1:length(nvals)
   plot(xvals,xvals.^nvals(nIdx),'LineWidth',2); 
end