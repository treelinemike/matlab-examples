% test central limit theorem sampling from exponential distribution
close all; clear all; clc;

lambda = 1/2;
n = 10000;

mu = 1/lambda;

xbars = [];

for sampleIdx = 1:100

s = exprnd(mu*ones(n,1));
xbar = mean(s);
sx = std(s);

xbars(end+1) = xbar;

end



sem_actual = std(xbars)
sem_expected = (1/lambda)/sqrt(n)
hist(xbars)
