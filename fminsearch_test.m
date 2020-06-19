% restart
close all; clear all; clc;

% define some observation data
data = [1 2 3];

% initial guess for optimal value of variable
x0 = 12;

% assemble observations and a variable for paramters
% into an anonymous function to be called by fminsearch
% note: necessary because... 
f = @(x)computeCost(data,x);

% run optimization and return result
[xOpt, finalCost] = fminsearch(f,x0)

% emperical evaluation of cost function
% x is the variable for which we are trying to compute an optimal value
% params are likely observed data
function err = computeCost(params,x)
    err = params(1)*x^2+params(2)*x+params(3);
    fprintf('evaluating at x=%+10.5f; cost = %+10.5f\n',x,err);
end