% restart
close all; clear all; clc;

% assemble observations and a variable for paramters
% into an anonymous function to be called by fminsearch
% note: necessary because... 
f = @(x)quad_cost(x);

% run optimization and return result
% least squares approach
options = optimset('Display','iter','PlotFcns',@optimplotfval);
[x_opt, cost_opt, exitflag, output] = fminsearch(f,12,options);


function cost = quad_cost(x)
cost = x^2;
end