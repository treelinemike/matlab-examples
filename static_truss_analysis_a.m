% Analyze system of equations for basic 7-member truss
% with one pinned and one rolling contact and a single applied load.
% This is the truss I'm using to review statics in ENGS 72 (Fall 2019,
% Winter 2020).
%
% Can double-check results with JHU applet: https://pages.jh.edu/~virtlab/bridge/bridge.htm

% restart
close all; clear; clc;

% truss angle = 60deg
theta = pi/3;

% coefficient/system matrix
A = [cos(theta) 0 0 0 0 1 0 1 0 0;
    sin(theta) 0 0 0 0 0 0 0 1 0;
    -cos(theta) 1 0 cos(theta) 0 0 0 0 0 0;
    -sin(theta) 0 0 -sin(theta) 0 0 0 0 0 0;
    0 -1 cos(theta) 0 -cos(theta) 0 0 0 0 0;
    0 0 -sin(theta) 0 -sin(theta) 0 0 0 0 0;
    0 0 -cos(theta) 0 0 0 -1 0 0 0;
    0 0 sin(theta) 0 0 0 0 0 0 1;
    0 0 0 -cos(theta) cos(theta) -1 1 0 0 0;
    0 0 0 sin(theta) sin(theta) 0 0 0 0 0];

% loading vector
b = zeros(10,1);
b(10) = 100;

% solve system
% note: negative values for member forces indicate compression
x = A\b