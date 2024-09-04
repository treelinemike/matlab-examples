% restart
close all; clear; clc;

R1 = 50; % [ohm]
R2 = 1; % [ohm]
C1 = 4.7e-6; % [F]
C2 = 1e-6; % [F]
A = [-1/(R2*C2) 1/(R2*C1); 1/(R2*C2) -1/(R2*C2)];
sysParams.A = A;

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 0.00001;      % [s] simulation end time
dt = 0.00000001;     % [s] timestep size

% initial conditions (state vector: [y y_dot]')
X0 = [4 2]'; % [m m/s]'
X = X0;

% data storage
time = t0;
data = X0;


vstar = -5:0.5:5;
vout = -6:0.5:6;
[X_VSTAR, Y_VOUT] = meshgrid(vstar,vout);
U = zeros(size(X_VSTAR));
V = zeros(size(X_VSTAR));

for voutIdx = 1:length(vout)
    for vstarIdx = 1:length(vstar)
        uv = A*[vstar(vstarIdx);vout(voutIdx)];
        %uv = uv/norm(uv);
        U(voutIdx,vstarIdx) = uv(1);
        V(voutIdx,vstarIdx) = uv(2);
    end
end

figure;
hold on; grid on; axis equal;
quiver(X_VSTAR,Y_VOUT,U,V,'Color',[0 0 0.8]);

[vec,val] = eig(A);
v1 = vec(:,1);
v2 = vec(:,2);
sf = (max(vstar)-min(vstar))/2;
plot(sf*[0 v1(1)], sf*[0 v1(2)],'-','LineWidth',2,'Color',[0.8 0 0]);
plot(sf*[0 v2(1)], sf*[0 v2(2)],'-','LineWidth',2,'Color',[0.8 0 0.8]);
xlabel('\bfv_{star}');
ylabel('\bfv_{out}');

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % apply bounce if needed
    if( X(1) < 0 && X(2) < 0 )
        X(2) = -1*sysParams.e*X(2);
    end
    
    % propagate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end
plot(data(1,:),data(2,:),'-','LineWidth',2,'Color',[0 0.8 0]);
plot(X0(1),X0(2),'.','MarkerSize',20,'Color',[0 0.8 0]);

figure;
ax(1) = subplot(2,1,1);
hold on; grid on;
plot(time,data(1,:),'-','LineWidth',1.6,'Color',[0 0.8 0]);
xlabel('\bfTime [s]');
ylabel('\bfv_{star}');

ax(2) = subplot(2,1,2);
hold on; grid on;
plot(time,data(2,:),'-','LineWidth',1.6,'Color',[0 0.8 0]);
xlabel('\bfTime [s]');
ylabel('\bfv_{out}');

linkaxes(ax,'x');

% propagate state
function Xdot = stateProp(t,X,sysParams)

% construct Xdot from differential equation
% note:     X    = [y y_dot]
% therefore Xdot = [y_dot y_ddot]
Xdot = sysParams.A*X;

end
