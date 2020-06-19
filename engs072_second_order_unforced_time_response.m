% restart
clear all; close all; clc;

% define mass, spring constant, damping coefficient
% CHANGE THESE PARAMETERS TO VARY THE SOLUTION
m = 5;
k = 2;
c = 2*m*sqrt(k/m);

% time vector
t = 0:0.01:100;

% initial conditions
y0 = 1;
ydot0 = 0;

% find roots of characteristic equation
s1 = -c/(2*m) + (1/(2*m))*sqrt(c^2-4*m*k);
s2 = -c/(2*m) - (1/(2*m))*sqrt(c^2-4*m*k);

% solve for linear combination coefficients
C = [1 1; s1 s2]\[y0; ydot0];

% assemble time response
y = C(1)*exp(s1*t) + C(2)*exp(s2*t);

% plot time response
figure;
hold on; grid on;
plot(t,y);
xlabel('\bfTime [s]');
ylabel('\bfy');