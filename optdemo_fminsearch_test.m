% restart
close all; clear all; clc;

% define parameters of cost function
% in this case coefficients a,b,c in a*x^2 + b*x + c
% could come from observed data
data = [1 2 3];

% initial guess for optimal value of variable
x0 = 4;

% assemble observations and a variable for paramters
% into an anonymous function to be called by fminsearch
% note: necessary because... 
f = @(x)computeCost(data,x);

% select optimization options
% allows performance plotting
opts = optimset('PlotFcn','optimplotfval');

% run optimization and return result
[xOpt, finalCost] = fminsearch(f,x0,opts);

% capture optimization performance figure
fh0 = gcf;
ax0 = gca;

% plot results
figure;
set(gcf,'Position',[0488 0361 7.362000e+02 2.376000e+02]);
ax1 = subplot(1,2,1);
hold on; grid on;
x = -5:0.01:5;
y = data(1)*x.^2+data(2)*x+data(3);
plot(x,y,'-','LineWidth',1.6,'Color',[0.0 0.0 0.8]);
y0 = computeCost(data,x0);
plot(x0,y0,'.','MarkerSize',30,'Color',[0 0.8 0]);
plot(xOpt,finalCost,'.','MarkerSize',30,'Color',[0.8 0 0]);
xlabel('\bfQuantity to Optimize');
ylabel('\bfCost');
legend('Cost','Initial Estimate','Optimal Point','Location','NorthWest');
ax2 = subplot(1,2,2);
copyobj(get(ax0,'children'),ax2);
hold on; grid on;
xlabel('\bfIteration');
ylabel('\bfCost');
close(fh0);
linkaxes([ax1,ax2],'y');

% emperical evaluation of cost function
% x is the variable for which we are trying to compute an optimal value
% params are likely observed data
function cost = computeCost(params,x)
    cost = params(1)*x^2+params(2)*x+params(3);
end