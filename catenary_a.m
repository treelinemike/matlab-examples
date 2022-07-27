% numerical optimization to find catenary
% based on shortest time between points in gravity with fixed energy

% restart
close all; clear; clc;

% parameters
xi = 5;
yi = 4;
xf = 15;
yf = 0;
N = 50;

% generate initial curve
x = linspace(xi,xf,N);
y_all = linspace(yi,yf,N);
y0 = y_all(2:end-1);

% create params struct
params.x = x;
params.g = 9.81;
params.yi = yi;
params.yf = yf;

% ensure initial height is above final height
assert(yi > yf,'Final height must be below initial height!');

% compute initial path time
t_initial = pathtime(y0,params)

% optimize path
opts = optimoptions('fmincon','StepTolerance',1e-16,'MaxFunctionEvaluations',1e6);
A = eye(length(y0));
b = params.yi*ones(length(y0),1);
[y_opt,t_min] = fmincon(@(y) pathtime(y,params),y0,A,b,[],[],[],[],[],opts);

% gather results
y = [params.yi y_opt params.yf];
t_min

% plot results
figure;
hold on; grid on;
plot(x,y_all,'LineWidth',1.6,'Color',[0.8 0 0]);
plot(x,y,'LineWidth',1.6,'Color',[0 0 0.8]);
daspect([1 1 1]);
legend(sprintf('Initial Path t = %0.2fs',t_initial),sprintf('Optimal Path t = %0.2fs',t_min));

% path time computation
function t = pathtime(y,params)
    y_all = [params.yi y params.yf];
    t = 0;
    v = zeros(size(y_all));
    for x_idx = 2:length(y_all)
        dy = y_all(x_idx)-y_all(x_idx-1);
        dx = params.x(x_idx)-params.x(x_idx-1);
        ds = sqrt((dx)^2+(dy)^2);        
        theta = atan2(dy,dx);
        a = -0.5*params.g*sin(theta);
        b = v(x_idx-1);
        c = -ds;
        these_roots = roots([a,b,c]);
        root_mask = these_roots > 0;
        this_t = min(these_roots(root_mask));
        v(x_idx) = sqrt(v(x_idx-1)^2-2*params.g*dy);
        assert(abs(v(x_idx) -  (v(x_idx-1)-params.g*sin(theta)*this_t)) < 1e-4,'Speed calc error!');
        t = t + this_t;
    end
end

