% NOTE: CONSIDER REMOVING LATEX INTERPRETER IN LEGENDS B/C MAC AND SOME
% OLDER VERSIONS OF MATLAB THROW ERRORS WITH THAT

% linear mass/spring/damper simulation

% restart
close all; clear all; clc;

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 4;         % [s] simulation end time
dt = 0.1;       % [s] timestep size

% initial conditions (state vector: [s s_dot]')
s_0     = 2;         % [m]
s_dot_0 = 0;         % [m/s]
X0 = [s_0 s_dot_0]'; % [m m/s]'
X = X0;

% system parameters
sysParams.m = 1;   % [kg]   mass
sysParams.k = 5;   % [N/m]  spring constant
sysParams.c = 3;   % [Ns/m] damping constant

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

% plot results
figure;
ah(1) = subplot(2,1,1);
hold on; grid on;
plot(time,zeros(size(time)),'k--','LineWidth',1.6);
plot(time,data(1,:),'b-','LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfPosition [m]');

ah(2) = subplot(2,1,2);
hold on; grid on;
plot(time,zeros(size(time)),'k--','LineWidth',1.6);
plot(time,data(2,:),'r-','LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfVelocity [m/s]');

linkaxes(ah,'y');

%% Animate result in a new plot
figure;
hold on; grid on;

% animate each frame of results
for tIdx = 1:size(data,2)
    
    % extract state at current timestep
    s = data(1,tIdx);
    s_dot = data(2,tIdx);
   
    % clear axes and start plotting the current frame
    cla;
    
    % plot XYZ (all black) and xyz (x=red, y=green, z=blue) coordinate frames
    plot(0,0,'k.','MarkerSize',50);
    plot(s,0,'k.','MarkerSize',50);
    plot([0 s],[0 0],'k-','LineWidth',6);
    plot(s,0,'o','MarkerSize',30,'LineWidth',6);%,'MarkerFaceColor',[1 1 1]);
    
    % finish formatting axes
    axis equal;
    xlabel('\bfX');
    ylabel('\bfY');
    xlim([-1 1]*1.2*max(abs(data(1,:))));
    ylim([-2 2]);
	drawnow;
end

% propagate state
function Xdot = stateProp(t,X,sysParams)

m = sysParams.m;
k = sysParams.k;
c = sysParams.c;

% deconstruct state vector
s     = X(1);
s_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [s s_dot]
% therefore Xdot = [s_dot s_ddot]
Xdot = zeros(2,1);
Xdot(1,:) = s_dot ;   % TODO: ENTER CORRECT EXPRESSION HERE
Xdot(2,:) = -(k/m)*s -(c/m)*s_dot ;   % TODO: ENTER CORRECT EXPRESSION HERE
end