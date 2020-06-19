% plot phase portraits for various parameter combinations of general second
% order, damped, unforced oscillator

function pendulum_phase_portrait_sim
% restart
close all; clear all; clc;
global sysParams;  % allow the system matrix to pass into the ODE call... not a good technique in general!
global doLinearize;
doLinearize = 1;

% define parameters for cases to analyze
sysParams.m = 1;  % pendulum acceleration is mass invariant
sysParams.l = 6;

% initial conditions of bottom plot [theta, theta_dot]'
ics = [5*pi/180 0]';

% time vector for trajectory simulation
tspan = [0 35];

% configure phase portrait grid
[x1,x2] = meshgrid(-6.5:0.5:6.5,-4:0.5:4);

% compute and display phase portrait


x1dot = x2;
if(doLinearize)
x2dot = -1*(9.81/sysParams.l)*x1;
else
    x2dot = -1*(9.81/sysParams.l)*sin(x1);
end

figure;
set(gcf,'Position',[4.882000e+02 2.034000e+02 5.744000e+02 5.584000e+02]);
subplot(8,1,1:5);
hold on; grid on;
quiver(x1,x2,x1dot,x2dot);
xlabel('\bfAngle [rad]');
ylabel('\bfAngular Velocity [rad/s]');
titleStr = 'Pendulum Phase Portrait';
if(doLinearize)
    titleStr = [titleStr ' - Linearized'];
end
title(titleStr);
axis equal;

% plot trajectory from ICs
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
[t,x] = ode45(@propDynamics,tspan,ics,opts);
x1traj = x(:,1)';
x2traj = x(:,2)';
plot(x1traj,x2traj,'-','Color',[0 0.7 0],'LineWidth',1.6);
set(gca,'XLim',[-8 8]);
set(gca,'YLim',[-5 5]);


% add fixed points
fixedPtLocs = [-2*pi 0; -pi 0; 0 0; pi 0; 2*pi 0]';
if(doLinearize)
    fixedPtLocs = zeros(size(fixedPtLocs));
end
plot(fixedPtLocs(1,1:2:size(fixedPtLocs,2)),fixedPtLocs(2,1:2:size(fixedPtLocs,2)),'ko','MarkerSize',5,'MarkerFace','k');
plot(fixedPtLocs(1,2:2:size(fixedPtLocs,2)),fixedPtLocs(2,2:2:size(fixedPtLocs,2)),'ko','MarkerSize',5);


% plot time series for detail trajectory
subplot(8,1,7:8);
hold on; grid on;
plot(t,0*ones(1,length(t)),'k--');
plot(t,x1traj,'Color',[0 0.7 0],'LineWidth',1.6);
xlabel('\bfTime');
ylabel('\bfx_1 (Position)');
set(gca,'YLim',[-4 4]);

% compute and display total energy in system
% to check accuracy of simulation
figure;
v = sysParams.l*x2traj;
E = 0.5*sysParams.m*v.^2 + sysParams.m*9.81*sysParams.l*(1-cos(x1traj));
plot(t,E);


% draw plot
drawnow;

end

function  xdot = propDynamics(t,x)
global sysParams;
global doLinearize;

theta = x(1);
theta_dot = x(2);
if(doLinearize)
    xdot = [theta_dot; -1*(9.81/sysParams.l)*theta];
else
    xdot = [theta_dot; -1*(9.81/sysParams.l)*sin(theta)];
end
end
