% NOTE: CONSIDER REMOVING LATEX INTERPRETER IN LEGENDS B/C MAC AND SOME
% OLDER VERSIONS OF MATLAB THROW ERRORS WITH THAT

% plot phase portraits for the simple pendulum

% restart
close all; clear all; clc;

% simulation time parameters
t0 = 0;       % [s] simulation start time
tf = 10;      % [s] simulation end time
dt = 0.1;     % [s] timestep size
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);

% define parameters for cases to analyze
sysParams = [];
sysParams.m = 1;  % pendulum acceleration is mass invariant
sysParams.l = 6;

% initial conditions (state vector: [theta theta_dot]')
theta_0     = 25*pi/180;      % [rad]
theta_dot_0 = 0;              % [rad/s]
X0 = [theta_0 theta_dot_0]';  % [rad rad/s]'
X = X0;

% data storage
time = [t0];
data = [X0];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % propigate state
    [T,X] = ode45(@(t,X) propDynamics(t,X,sysParams),odeTime,X,opts);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% configure phase portrait grid
[x1,x2] = meshgrid(-6.5:0.5:6.5,-4:0.5:4);

% compute phase portrait
x1dot = x2;
x2dot = -1*(9.81/sysParams.l)*sin(x1);

% plot phase portrait
figure;
set(gcf,'Position',[4.882000e+02 2.034000e+02 5.744000e+02 5.584000e+02]);
subplot(8,1,1:5);
hold on; grid on;
quiver(x1,x2,x1dot,x2dot);
xlabel('\bfAngle [rad]');
ylabel('\bfAngular Velocity [rad/s]');
title('Pendulum Phase Portrait');
axis equal;

% add fixed points
fixedPtLocs = [-2*pi 0; -pi 0; 0 0; pi 0; 2*pi 0]';
plot(fixedPtLocs(1,1:2:size(fixedPtLocs,2)),fixedPtLocs(2,1:2:size(fixedPtLocs,2)),'ko','MarkerSize',5,'MarkerFace','k');
plot(fixedPtLocs(1,2:2:size(fixedPtLocs,2)),fixedPtLocs(2,2:2:size(fixedPtLocs,2)),'ko','MarkerSize',5);

% plot trajectory from simulation
x1traj = data(1,:)';
x2traj = data(2,:)';
plot(x1traj(1),x2traj(1),'.','Color',[0 0.7 0],'MarkerSize',20);
plot(x1traj,x2traj,'-','Color',[0 0.7 0],'LineWidth',1.6);
set(gca,'XLim',[-8 8]);
set(gca,'YLim',[-5 5]);

% plot time series for detail trajectory
subplot(8,1,7:8);
hold on; grid on;
plot(time,0*ones(1,length(time)),'k--');
plot(time,x1traj,'Color',[0 0.7 0],'LineWidth',1.6);
xlabel('\bfTime');
ylabel('\bfx_1 (Position)');
set(gca,'YLim',[-4 4]);


%% Animate result in a new plot
figure;
hold on; grid on;

% animate each frame of results
for tIdx = 1:size(data,2)
    
    % extract state at current timestep
    theta = data(1,tIdx);
    theta_dot = data(2,tIdx);
   
    % recover length
    l = sysParams.l;
    
    % determine bob location
    bob_x = l*sin(theta);
    bob_y = -l*cos(theta);
    
    % clear axes and start plotting the current frame
    cla;
    
    % plot XYZ (all black) and xyz (x=red, y=green, z=blue) coordinate frames
    plot(0,0,'k.','MarkerSize',30);
    plot([0 bob_x],[0 bob_y],'k-','LineWidth',6);
    plot(bob_x,bob_y,'o','MarkerSize',30,'LineWidth',6,'MarkerFaceColor',[1 1 1]);
    
    % finish formatting axes
    axis equal;
    xlabel('\bfX');
    ylabel('\bfY');
    xlim([-1.2*l 1.2*l]);
    ylim([-1.2*l 1.2*l]);
	drawnow;
end

% propagate state
function  Xdot = propDynamics(t,X,sysParams)

% recover paramters
m = 1;  % pendulum acceleration is mass invariant
l = 6;

% deconstruct state vector
theta = X(1);
theta_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [theta      theta_dot]
% therefore Xdot = [theta_dot  theta_ddot]
Xdot = zeros(2,1);
Xdot(1,:) = theta_dot;
Xdot(2,:) = -1*(9.81/l)*sin(theta);

end
