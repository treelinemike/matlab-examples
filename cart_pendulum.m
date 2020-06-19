% cart/pendulum simulation
% equations of motion from Newton or Lagrange
function cart_pendulum
% restart
close all; clear all; clc;

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 4;         % [s] simulation end time
dt = 0.015;     % [s] timestep size
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);

% initial conditions (state vector: [x x_dot theta theta_dot]')
x_0         = 0;  % [m]
x_dot_0     = 0;  % [m/s]
theta_0     = 160*pi/180;  % [rad]
theta_dot_0 = 0;  % [rad/s]
X0 = [x_0 x_dot_0 theta_0 theta_dot_0]'; % [m m/s rad rad/s]'
X = X0;

% system parameters
sysParams.ma = 2;    % [kg] cart mass
sysParams.mb = 3;    % [kg] bob mass
sysParams.l  = 0.7;  % [m]  pendulum length
sysParams.g = 9.81;  % [m/s^2]  acceleration of gravity

% data storage
time = [t0];
data = [X0];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % propigate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X,opts);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% compute total energy in system
xa = data(1,:);
x_dot = data(2,:);
theta = data(3,:);
theta_dot = data(4,:);
ma = sysParams.ma;
mb = sysParams.mb;
l = sysParams.l;
g = sysParams.g;
T = 0.5*ma*x_dot.^2 + 0.5*mb*(x_dot + l*theta_dot.*cos(theta)).^2 + 0.5*mb*(l*theta_dot.*sin(theta)).^2;
V = mb*g*l*(1-cos(theta));
E = T+V;

% plot total energy
figure;
hold on; grid on;
plot(time,E-E(1),'-','LineWidth',1.6,'Color',[0 0.7 0]);
xlabel('\bfTime [sec]');
ylabel('\bf\Delta Energy [J]');
title('\bfEnergy Check','FontSize',12);

%% animate result in a new plot

% initialize figure
figure;
hold on; grid on;

% define cart patch geometry
cart.v = [ -0.15 0;
    -0.15 0.15;
    0.15 0.15;
    0.15 0];
cart.f = [1 2 3 4 1];

% animate each frame of results
for tIdx = 1:size(data,2)
    
    % extract state at current timestep
    y_shift = 0.075;
    xa = data(1,tIdx);
    x_dot = data(2,tIdx);
    theta = data(3,tIdx);
    theta_dot = data(4,tIdx);
    ya = y_shift;
    xb = xa+sysParams.l*sin(theta);
    yb = ya-sysParams.l*cos(theta);
    cart.thisv = cart.v + xa*repmat([1,0],4,1);
    
    % clear axes and start plotting the current frame
    cla;
    
    % plot system
    plot([-3 3],[0 0],'k-','LineWidth',2);
    patch('Faces',cart.f,'Vertices',cart.thisv,'FaceColor','flat','EdgeColor','k','LineWidth',2,'FaceColor',0.7*ones(1,3));
    plot(xa,ya,'k.','MarkerSize',50);
    plot([xa xb],[ya yb],'k-','LineWidth',4);%,'MarkerFaceColor',[1 1 1]);
    plot(xb,yb,'b.','MarkerSize',100);%,'MarkerFaceColor',[1 1 1]);
    
    % compute and display center of mass
    xcm = (sysParams.ma*xa + sysParams.mb*xb)/(sysParams.ma+sysParams.mb);
    plot(xcm,0,'r.','MarkerSize',20);
    
    % finish formatting axes
    axis equal;
    xlabel('\bfX');
    ylabel('\bfY');
    xlim([-1 1]);
    ylim([-0.8 1]);
    title('\bfPendulum/Cart Simulation','FontSize',12);
    drawnow;
end
end

% propagate state
function Xdot = stateProp(t,X,sysParams)

% recover parameters
ma = sysParams.ma;
mb = sysParams.mb;
l = sysParams.l;
g = sysParams.g;

% deconstruct state vector
x         = X(1);
x_dot     = X(2);
theta     = X(3);
theta_dot = X(4);

% solve for x_ddot and theta_ddot

% LAGRANGE u = [x_ddot, theta_ddot]'
A = [ ma+mb,          mb*l*cos(theta);
    l*cos(theta),   l^2  ];
b = [ mb*l*theta_dot^2*sin(theta);
    -g*l*sin(theta)];
u = A\b;

% % NEWTON  u = [x_ddot, theta_ddot, Fn]'
% A = [ ma,             0,                -sin(theta);
%       0,              mb*l*sin(theta),  -cos(theta);
%       l*cos(theta),   l^2,              0];
% b = [ 0;
%      -mb*g - mb*l*(theta_dot^2)*cos(theta);
%      -g*l*sin(theta)];
% u = A\b;

% construct Xdot from differential equation
% note:     X    = [x x_dot theta theta_dot]'
% therefore Xdot = [x_dot x_ddot theta_dot theta_ddot]'
Xdot = zeros(4,1);
Xdot(1,:) = x_dot;
Xdot(2,:) = u(1);
Xdot(3,:) = theta_dot;
Xdot(4,:) = u(2);
end