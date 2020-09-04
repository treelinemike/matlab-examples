% Based on ENGS 199-01 HW #4 - Problem #1
% Generate observed ata 

% restart
close all; clear; clc;

grad_type = 'stochastic';   %'descent','coord','stochastic'

% parameters
numobs = 20;
mu = [2,3];
varx = 1;
vary = 3;
mLim = [-3 7];
bLim = [-3 8];
mStep = 0.5;
bStep = 0.5;
gStep = 0.001;
plotIters = [1, 2, 5, 10, 100];

% generate observations
covxy_max = sqrt(varx*vary);   % corresponds to corr(x,y) = 1; limit here is sqrt(3) = 1.7321
% covxy = covxy_max;
covxy = 1.5;
sigma = [varx,covxy;covxy,vary];
rng default;
r = mvnrnd(mu,sigma,numobs);
obs_x = r(:,1); obs_y = r(:,2);

% compute least squares solution
H = [ones(size(obs_x)) obs_x];
xLS = inv(H'*H)*H'*obs_y;   % alternatively xLS = H\obs_y
JLS = (obs_y-H*xLS)'*(obs_y-H*xLS);

% compute error surface (i.e. cost function)
m = mLim(1):mStep:mLim(end);
b = bLim(1):bStep:bLim(end);
J = ones(length(b),length(m));
for thisBIdx = 1:length(b)
    for thisMIdx = 1:length(m)
        thisB = b(thisBIdx);
        thisM = m(thisMIdx);
        thisX = [thisB thisM]'; 
        J(thisBIdx,thisMIdx) = (obs_y-H*thisX)'*(obs_y-H*thisX);
    end
end
[Jm,Jb] = gradient(J); % used to construct quiver plot

% plot observations
for i = 1:2
    figure;
    hold on; grid on;
    ph0 = plot(obs_x,obs_y,'r.','MarkerSize',20);
    xlabel('\bfObservation X Component');
    ylabel('\bfObservation Y Component');
    axis equal;
end

% plot least squares solution
yLS = H*xLS;
ph0(end+1) = plot(obs_x,yLS,'--','Color',[0 0.7 0],'LineWidth',1.6);
legStr = {'Observations','Least Squares Solution'};

% initial conditions
x0 = [-1 5]';
x0 = [6 6]';
J0 = (obs_y-H*x0)'*(obs_y-H*x0);
data = [x0' J0];

% now we repurpose x to be our vector of unknown states [b, m]
x = x0;

% iterate through gradient descent in a rudimentary fashion
% note: this implementation does NOT change the step size / learning rate
% (gStep) but ideally it would do a line/ray search to optimize the step size
% at each iteration
for iter = 1:5000
    
    % compute the local gradient
    % if we didn't have an analytical expression
    % we could do this emperically by takings steps
    % and approximating the local gradient
    grad = -2*H'*obs_y+2*H'*H*x; % [b,m]' ... just the gradient of the cost function J = (Hx-y)'*(Hx-y)
    
    % change gradient if using coordinate or stochastic gradient descent
    switch grad_type
        case 'descent'
            grad_adj = grad;
        case 'coordinate'
            grad_adj = grad.*[mod(iter,2) mod(iter+1,2)]';
        case 'stochastic'
            theta = pi*rand;
            v = [cos(theta); sin(theta)];
            grad_adj = dot(grad,v)*v;
        otherwise
            error('Invalid gradient type!');
    end
    
    % take a step and compute the cost at the new location
    x = x - gStep * grad_adj;
    thisJ = (obs_y-H*x)'*(obs_y-H*x);
    
    % store data
    data(end+1,:) = [x' thisJ]; % [b,m,cost]
    
    % plot estimated line at predefined iterations
    if( max(iter == plotIters) )
        thisY = H*x;
        ph0(end+1) = plot(obs_x,thisY,'LineWidth',1.6);
        legStr{end+1} = sprintf('Fit at Step #%d',iter);
    end
end

% add legend to plot
legend(ph0,legStr,'Location','NorthWest');
set(gca,'XLim',[-8 10]);
set(gca,'YLim',[-3 12]);

% plot gradient descent surface and path
figure;
set(gcf,'Position',[0481 6.340000e+01 5.032000e+02 7.096000e+02]);
t = tiledlayout(2,1);
t.Padding = 'none';  % 'normal', 'compact', or 'none'
t.TileSpacing = 'none';  % 'normal', 'compact', or 'none'
nexttile(1);
hold on; grid on;
sh = surf(m,b,J);
contour(m,b,J,40);
quiver(m,b,Jm,Jb);
ph1 = plot3(x0(2),x0(1),J0,'.','MarkerSize',25,'Color',[0 0.9 0]);
plot3(x0(2),x0(1),J0,'o','MarkerSize',10,'Color',[0 0 0],'LineWidth',2);
ph1(end+1) = plot3(xLS(2),xLS(1),JLS,'.','MarkerSize',25,'Color',[0.8 0.0 0.8]);
plot3(xLS(2),xLS(1),JLS,'o','MarkerSize',10,'Color',[0 0 0],'LineWidth',2);
ph1(end+1) = plot3(data(:,2),data(:,1),data(:,3),'r.-','LineWidth',3);
xlabel('\bfSlope');
ylabel('\bfIntercept');
zlabel('\bfCost');
% legend(ph1,{'Initial Condition','Least Squares Solution','Gradient Descent'});
view([-127,36]);
set(sh,'FaceAlpha',0.6);
set(gca,'DataAspectRatio',[1 1 1500])

% show same data with contour plot only
nexttile(2);
hold on; grid on;
contour(b,m,J',40);
quiver(b,m,Jb',Jm');
plot(x0(1),x0(2),'.','MarkerSize',25,'Color',[0 0.9 0]);
plot(x0(1),x0(2),'o','MarkerSize',10,'Color',[0 0 0],'LineWidth',2);
plot(xLS(1),xLS(2),'.','MarkerSize',25,'Color',[0.8 0.0 0.8]);
plot(xLS(1),xLS(2),'o','MarkerSize',10,'Color',[0 0 0],'LineWidth',2);
plot(data(:,1),data(:,2),'r.-','LineWidth',3);
xlabel('\bfIntercept');
ylabel('\bfSlope');
set(gca,'XDir','reverse');
axis equal;

% plot cost function vs. iteration
figure;
set(gcf,'Position',[ 610   428   710   335]);
loglog((1:size(data,1))',data(:,3),'r-','LineWidth',1.6);
grid on;
xlabel('\bfIteration');
ylabel('\bfCost Function Value');

% save cost trajectory
costTraj = [(1:size(data,1))' data(:,3)];
save(['costTraj' grad_type '.mat'],'costTraj');


%% plot cost function vs. iteration
figure;
set(gcf,'Position',[0610 4.722000e+02 4.718000e+02 2.908000e+02]);

load('costTrajdescent.mat');
ph = loglog(costTraj(:,1),costTraj(:,2),'-','LineWidth',1.6,'Color',[0.8 0.0 0.0]);
hold on; grid on;

load('costTrajcoordinate.mat');
ph(end+1) = loglog(costTraj(:,1),costTraj(:,2),'-','LineWidth',1.6,'Color',[0.0 0.8 0.0]);

load('costTrajstochastic.mat');
ph(end+1) = loglog(costTraj(:,1),costTraj(:,2),'-','LineWidth',1.6,'Color',[0.0 0.0 0.8]);

xlabel('\bfIteration');
ylabel('\bfCost');
legend(ph,{'Gradient','Coordinate','Stochastic'});

