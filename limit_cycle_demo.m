% LIMIT CYCLE PHASE PORTRAIT EXAMPLE
% Quiver plotting fixed by Claude

% restart
close all; clear all; clc;

% general options
params.alpha = 0.28;   % critical value ~= 0.19
params.Rp = 50;    % [Ohm]
params.C = 1e-4;  % [F]
params.L = 0.1;  % [H]
params.cubic_scale = 0.1;

% initial conditions for trajectory
VC0 = 0.45; % [V]
IL0 = 0;    % [A]

% simulation time parameters
t0 = 0;        % [s] simulation start time
tf = 1;        % [s] simulation end time
dt = 0.0001;   % [s] timestep size
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

% compute phase portrait
% x1 = V_C
% x2 = I_L
[x1,x2] = meshgrid(-0.5:0.05:0.5,-0.02:0.001:0.02);
x1dot = -(1/params.C)*x2 - (1/(params.Rp*params.C))*x1 + params.cubic_scale*(- (1/params.C)*x1.^3 + (1/params.C)*params.alpha*x1);
x2dot = (1/params.L)*x1;

% ---- fix BOTH direction and arrowhead skew ----
% This is a Claude contribution...
% Problem: quiver draws arrowhead geometry in *data units*. With V spanning
% 0.2 and A spanning 0.006 (~33x mismatch), the renderer stretches the
% arrowhead triangle by that same ratio -> skewed/lopsided heads, even
% after fixing arrow direction/length alone.
%
% Fix: rescale the current axis internally so its numeric span matches the
% voltage axis span, plot everything in that space, then relabel the
% y-ticks back to true Amps. Displayed axis values are unaffected -- only
% the internal plotting coordinate is rescaled.

xspan = range(x1(:));   % V range
yspan = range(x2(:));   % A range
s = xspan / yspan;      % scale factor: 1 "unit" of A now spans same numeric range as V

x2s    = x2 * s;         % rescaled current, for plotting only
x2dots = x2dot * s;      % rescaled current derivative, for plotting only

% normalize direction/length in this now-isotropic (equal-span) space
mag = sqrt(x1dot.^2 + x2dots.^2);
mag(mag == 0) = 1; % avoid divide-by-zero at fixed point
x1dot_n  = x1dot  ./ mag;
x2dots_n = x2dots ./ mag;

arrowFrac = 0.04;   % arrow length as fraction of axis span; tune to taste
x1dot_plot = x1dot_n  * xspan * arrowFrac;
x2dot_plot = x2dots_n * xspan * arrowFrac;  % xspan on purpose: space is isotropic now

% plot phase portrait
figure(1);
quiver(x1,x2s,x1dot_plot,x2dot_plot,0);  % scale=0: use our own pre-scaled lengths
hold on; grid on;
xlabel('\bfCapacitor Voltage [V]');
ylabel('\bfInductor Current [A]');
title('Phase Portrait');

% relabel y-ticks to show true current values (divide back out of scaled space)
yt = yticks;
yticklabels(compose('%.4f', yt / s));

% LSIM()
% A = [-1/(params.Rp*params.C) -1/params.C; 1/params.L 0];
% B = [0; 0];
% C = eye(2);
% D = [0; 0];
% sys = ss(A,B,C,D);
% t = 0:0.0001:0.1;
% [y,t_out] = lsim(sys,[],t,[0.1;0]);

% data storage
X0 = [VC0 IL0]';
X = X0;
time = [t0];
data = [X0];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];

    % propigate state
    [T,X] = ode45(@(t,X) propDynamics(t,X,params),odeTime,X,opts);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% add simulated trajectory to phase portrait
figure(1);
plot(data(1,:),data(2,:)*s,'LineWidth',1.6,'Color',[0.8 0 0]);  % *s to match rescaled quiver space

% propagate state
function  Xdot = propDynamics(t,X,params)

% recover paramters
C = params.C;  % note: pendulum acceleration is mass invariant
L = params.L;
Rp = params.Rp;
alpha = params.alpha;
cubic_scale = params.cubic_scale;

% deconstruct state vector
VC = X(1);
IL = X(2);

% construct Xdot from differential equation
Xdot = zeros(2,1);
Xdot(1,:) = -(1/C)*IL - (1/(Rp*C))*VC + cubic_scale*(-(1/C)*VC.^3 + (1/C)*alpha*VC);
Xdot(2,:) = (1/L)*VC;

end