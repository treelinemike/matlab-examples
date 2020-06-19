% restart
close all; clear all; clc;

% data storage
leg_plots = [];
leg_strings = {};

dx = 0.01;

% create actual function that we'll try to approximate
x1 = 0:dx:2*pi-pi/2;
y1 = cumtrapz(sin(x1)+1);
x2 = x1(end)+dx:dx:4*pi-pi/2;
y2 = cumtrapz(2*(sin(x2)+1)) + y1(end);
x3 = x2(end)+dx:dx:16;
y3 = y2(end)*ones(1,length(x3));
x = [x1 x2 x3]';
y = [y1 y2 y3]';

figure;
% set(gcf,'Position',[1.686600e+03 1.482000e+02 1.102400e+03 0388]);
hold on; grid on;
leg_plots(end+1) = plot(x,y,'LineWidth',6,'Color',[0 1 1]);
leg_strings{end+1} = 'Sample Data';
title('\bfLinear Function Approximation');

% Linear Regression
A = [ones(size(x,1),1), x];
alpha_hat = (A'*A)\(A'*y);
y_hat = A*alpha_hat;
leg_plots(end+1) = plot(x,y_hat,'r','LineWidth',1.6);
leg_strings{end+1} = 'Standard Linear Model';
legend(leg_plots,leg_strings);

% Fourier Series Approximation
n_max = 20;
n_vals = -n_max:1:n_max;
i = sqrt(-1);
L = max(x) - min(x);
A = exp(i*2*pi.*n_vals.*x/L);
alpha_hat = (A'*A)\(A'*y);
y_hat = A*alpha_hat;
leg_plots(end+1) = plot(x,y_hat,'Color',[0 0.7 0],'LineWidth',1.6);
leg_strings{end+1} = 'Fourier Basis';
legend(leg_plots,leg_strings,'Location','SouthEast');

% Radial Basis Functions
n = 20;
sigma = 2;
mu = min(x): (max(x)-min(x))/(n-1) : max(x);
A = exp(-(x-mu).^2/sigma);
alpha_hat = (A'*A)\(A'*y);
y_hat = A*alpha_hat;
leg_plots(end+1) = plot(x,y_hat,'Color',[1 0 1],'LineWidth',1.6);
leg_strings{end+1} = 'Radial Basis Functions';
legend(leg_plots,leg_strings,'Location','SouthEast');


