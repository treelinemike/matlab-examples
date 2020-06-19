% plot phase portraits for various parameter combinations of general second
% order, damped, unforced oscillator
% written by m.kokko, 11-jan-17

function second_order_unforced

% restart
close all; clear all; clc;
global A;

% options 
doSaveFrames = 0;

% define parameters for cases to analyze
paramVals = [0 1 1; 0.25 1 1; 0.5 1 1; 1 1 1; 1.9 1 1; 2 1 1; 2.1 1 1; 10 1 1; 50 1 1];
% paramVals = [ 0 logspace(-1,2,60) ]
% paramVals = [paramVals; ones(2,length(paramVals))]';

% initial conditions
dtheta = pi/12;
theta = 0:dtheta:2*pi-dtheta;
ics = [cos(theta); sin(theta)];
ics = [ics [3 3; -3 3; 3 -3; -3 -3]'];
detail_ic = [0 3]';

% time vector for trajectory simulation
tsim = 0:0.1:35;

% configure phase portrait grid
[x1,x2] = meshgrid(-6:1:6,-4:1:4);

% compute and display phase portrait for each set of parameters
for idx = 1:size(paramVals,1)
    
    b = paramVals(idx,1);
    k = paramVals(idx,2);
    m = paramVals(idx,3);
    A = [ 0 1; -k/m -b/m ];
    
    % classify
    if( b == 0 )
        classStr = 'Undamped';
    elseif( b^2 < 4*m*k )
        classStr = 'Underdamped';
    elseif( b^2 == 4*m*k )
        classStr = 'Critically Damped';
    elseif( b^2 > 4*m*k )
        classStr = 'Overdamped';
    end
    
    [vec,val] = eig(A);
    x1dot = A(1,1).*x1+A(1,2).*x2;
    x2dot = A(2,1).*x1+A(2,2).*x2;
    
    figure;
    subplot(8,1,1:5);
    hold on; grid on;
    quiver(x1,x2,x1dot,x2dot);
    xlabel('\bfx_1 (Position)');
    ylabel('\bfx_2 (Velocity)');
    title(sprintf([classStr ': b=%0.2f, k=%0.2f, m=%0.2f'],b,k,m));
    axis equal;
    
    % plot trajectories from ICs
    for icidx = 1:size(ics,2)
        [~,x] = ode45(@propDynamics,tsim,ics(:,icidx));
        x1traj = x(:,1)';
        x2traj = x(:,2)';
        plot(x1traj,x2traj,'r');

    end
    
    % plot last trajectory for which time series trace will be shown
    [t,x] = ode45(@propDynamics,tsim,detail_ic);
    x1traj = x(:,1)';
    x2traj = x(:,2)';
    plot(x1traj,x2traj,'-','Color',[0 0.7 0],'LineWidth',1.6);
    set(gca,'XLim',[-8 8]);
    set(gca,'YLim',[-5 5]);
    
    % draw eigenvectors
    if(isreal(vec(:,1)))
        plot([ 0 vec(1,1)],[0 vec(2,1)],'m','LineWidth',1.6);
    else
        plot([ 0 real(vec(1,1))],[0 real(vec(2,1))],'c','LineWidth',1.6);
    end
    
    if(isreal(vec(:,1)))
        plot([ 0 vec(1,2)],[0 vec(2,2)],'m','LineWidth',1.6);
    else
        plot([ 0 real(vec(1,2))],[0 real(vec(2,2))],'c','LineWidth',1.6);
    end
    
    % add fixed points
    fixedPtLocs = [0 0]';
    fixedPtLocs = [fixedPtLocs null(A)];
    plot(fixedPtLocs(1,:),fixedPtLocs(2,:),'k.','MarkerSize',20);
    
    % plot time series for detail trajectory
    subplot(8,1,7:8);
    hold on; grid on;
    plot(t,0*ones(1,length(t)),'k--');
    plot(t,x1traj,'Color',[0 0.7 0],'LineWidth',1.6);
    xlabel('\bfTime');
    ylabel('\bfx_1 (Position)');
    set(gca,'YLim',[-4 4]);
    
    % draw plot
    drawnow;
    
    % save figure;
    if(doSaveFrames)
        saveas(gcf,sprintf('frame%03d.png',idx));
        close all;
    end
end
end

function  xdot = propDynamics(t,x)
global A;

% LTI so we don't use t
xdot = A*x;

end
