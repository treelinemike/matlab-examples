% gravity simulation

% restart
close all; clear all; clc;

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 4;         % [s] simulation end time
dt = 0.05;     % [s] timestep size

% initial conditions (state vector: [y y_dot]')
y_0     = 0;         % [m]
y_dot_0 = 0;         % [m/s]
X0 = [y_0 y_dot_0]'; % [m m/s]'
X = X0;

% system parameters
sysParams.m = 1;    % [kg]     mass
sysParams.g = 9.81; % [m/s^2]  acceleration of gravity

% data storage
time = [t0];
data = [X0];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % propigate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% compute total energy in system
y = data(1,:);
y_dot = data(2,:);
E = 0.5*sysParams.m*y_dot.^2+sysParams.m*sysParams.g*y;

% plot results
figure;
ah(1) = subplot(3,1,1);
hold on; grid on;
plot(time,zeros(size(time)),'k--','LineWidth',1.6);
plot(time,data(1,:),'b-','LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfPosition [m]');

ah(2) = subplot(3,1,2);
hold on; grid on;
plot(time,zeros(size(time)),'k--','LineWidth',1.6);
plot(time,data(2,:),'r-','LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfVelocity [m/s]');

ah(3) = subplot(3,1,3);
hold on; grid on;
plot(time,E,'-','LineWidth',1.6,'Color',[0 0.7 0]);
xlabel('\bfTime [sec]');
ylabel('\bfEnergy [J]');

linkaxes(ah,'x');

%% Animate result in a new plot
figure;
hold on; grid on;

% animate each frame of results
for tIdx = 1:size(data,2)
    
    % extract state at current timestep
    y = data(1,tIdx);
    y_dot = data(2,tIdx);
   
    % clear axes and start plotting the current frame
    cla;
    
    % plot XYZ (all black) and xyz (x=red, y=green, z=blue) coordinate frames
    plot(0,0,'k.','MarkerSize',50);
    plot(0,y,'k.','MarkerSize',50);
%     plot([0 0],[0 y],'k-','LineWidth',6);
    plot(0,y,'o','MarkerSize',30,'LineWidth',6);%,'MarkerFaceColor',[1 1 1]);
    
    % finish formatting axes
    axis equal;
    xlabel('\bfX');
    ylabel('\bfY');
    xlim([-20 20]);
    ylim([-1 0]*1.2*max(abs(data(1,:))));
	drawnow;
end

% propagate state
function Xdot = stateProp(t,X,sysParams)

% recover parameters
m = sysParams.m;
g = sysParams.g;

% deconstruct state vector
y     = X(1);
y_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [y y_dot]
% therefore Xdot = [y_dot y_ddot]
Xdot = zeros(2,1);
Xdot(1,:) = y_dot;
Xdot(2,:) = -1*g;  
end