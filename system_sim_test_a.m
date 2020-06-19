% Show that differential equation can be framed in terms of either physical
% state variables (position, velocity), canonical state variables
% (system decoupled along its eigenvectors -> x_1 and x_2), or a linear combination of
% canonical state variables that results in a real system (x_3 and x_4).
%
% NOTE: To run the decoupled simulation the new variables need to be
% broken up into real and imaginary components; for some reason ODE45
% doesn't appear to work properly with the complex system otherwise.
%
% Simple unforced, undamped linear mass/spring system.

% restart
close all; clear all; clc;

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 4;         % [s] simulation end time
dt = 0.1;       % [s] timestep size

% system parameters
sysParams.m = 1;        % [kg]   mass
sysParams.k = 9.8696;   % [N/m]  spring constant

% state space representation using PHYSICAL state vars x = [pos vel]
sysParams.A = [0 1; -sysParams.k/sysParams.m 0];   % system matrix
B = 0;               % no forcing function

% initial conditions (state vector: [s s_dot]')
s_0     = 2;         % [m]
s_dot_0 = 0;         % [m/s]
X0 = [s_0 s_dot_0]'; % [m m/s]'
X = X0;

% state space representation using CANONICAL state vars (eigs)
[S,L] = eig(sysParams.A);
sysParams.A_star = L; % == inv(S)*A*S;
X_star0 = S\X0;
X_star = X_star0;

% state space representation using transformation of canonical state vars
% s.t. the system is real; eigenvectors combined s.t. system becomes real
% in a system with other real eigenvalues this would be the only thing that
% would need to change
T = [0.5 -0.5*i; 0.5 0.5*i]; % from Maybeck: Stochastic Models Estimation and Control vol 1 pg. 34
sysParams.A_dstar = inv(T)*L*T;
X_dstar0 = T\(S\X0);
X_dstar = X_dstar0;


% set initial states
X = [X;X_star;X_dstar];
X = [X(1:2); real(X(3)); imag(X(3)); real(X(4)); imag(X(4)); X(5:6)];

% data storage
time = [t0];
data = [X];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % propigate state
    [Tout,Xout] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X);
    X = Xout(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = Tout(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% 
X_star_expected = inv(S)*data(1:2,:);
X_star = [data(3,:)+i*data(4,:); data(5,:)+i*data(6,:)];
X_recon = S*X_star;
X_dstar_expected = inv(T)*inv(S)*data(1:2,:);
X_dstar = data(7:8,:);
X_dstar_recon = S*T*X_dstar;

% plot results
figure;
ah = subplot(6,1,1);
hold on; grid on;
plot(time,real(X_recon(1,:)),'b-','LineWidth',1.6);
plot(time,real(X_dstar_recon(1,:)),'m:','LineWidth',2.5);
plot(time,data(1,:),'k--','LineWidth',1.6);
legend('\bfFrom x_1 & x_2','\bfFrom x_3 & x_4','\bfTrue');
xlabel('\bfTime [sec]');
ylabel('\bfPosition [m]');

ah(end+1) = subplot(6,1,2);
hold on; grid on;
plot(time,real(X_recon(2,:)),'b-','LineWidth',1.6);
plot(time,real(X_dstar_recon(2,:)),'m:','LineWidth',2.5);
plot(time,data(2,:),'k--','LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfVelocity [m/s]');

ah(end+1) = subplot(6,1,3);
hold on; grid on;
plot(time,real(X_star(1,:)),'b-','LineWidth',1.6);
plot(time,imag(X_star(1,:)),'r-','LineWidth',1.6);
plot(time,real(X_star_expected(1,:)),'k--','LineWidth',2);
plot(time,imag(X_star_expected(1,:)),'k:','LineWidth',2);
legend('\bfRe','\bfIm','\bfRe (true)','\bfIm (true)');
xlabel('\bfTime [sec]');
ylabel('\bfx_1');

ah(end+1) = subplot(6,1,4);
hold on; grid on;
plot(time,real(X_star(2,:)),'b-','LineWidth',1.6);
plot(time,imag(X_star(2,:)),'r-','LineWidth',1.6);
plot(time,real(X_star_expected(2,:)),'k--','LineWidth',2);
plot(time,imag(X_star_expected(2,:)),'k:','LineWidth',2);
xlabel('\bfTime [sec]');
ylabel('\bfx_2');

ah(end+1) = subplot(6,1,5);
hold on; grid on;
plot(time,X_dstar(1,:),'b-','LineWidth',1.6);
plot(time,real(X_dstar_expected(1,:)),'k--','LineWidth',2);
xlabel('\bfTime [sec]');
ylabel('\bfx_3');

ah(end+1) = subplot(6,1,6);
hold on; grid on;
plot(time,X_dstar(2,:),'b-','LineWidth',1.6);
plot(time,real(X_dstar_expected(2,:)),'k--','LineWidth',2);
xlabel('\bfTime [sec]');
ylabel('\bfx_4');


linkaxes(ah,'y');


% propagate state
function Xdot = stateProp(t,X,sysParams)
    thisX = [X(1:2); X(3)+i*X(4); X(5)+i*X(6); X(7:8)];
    thisXdot = blkdiag(sysParams.A,sysParams.A_star,sysParams.A_dstar)*thisX;
    Xdot = [thisXdot(1:2); real(thisXdot(3)); imag(thisXdot(3)); real(thisXdot(4)); imag(thisXdot(4)); thisXdot(5:6)];
end