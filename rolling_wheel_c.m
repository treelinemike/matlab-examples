% restart
close all; clear all; clc;

% wheel radius
r0 = 1;

% maximum trajectory x position
trajMaxX = 3.2;

% constant angular velocity
theta_dot = -1; % [rad/s]
theta_ddot = 0; % [rad/s^2]

% time
t=0:pi/16:2*pi;
t_samp = 0:pi/4:2*pi;

% generate angular positions
theta = theta_dot*t;
theta_samp = theta_dot*t_samp;

% initialize plot
figure;
set(gcf,'Position',[3.602000e+02 4.180000e+01 8.056000e+02 7.408000e+02]);
hold on;

% draw circle
circTheta = 0:0.01:2*pi;
circX = r0*cos(circTheta);
circY = r0+r0*sin(circTheta);
plot(circX,circY,'b-','LineWidth',1.6);
plot(0,r0,'k.','MarkerSize',15);

vscale = 2;
ascale = 1.5;
axHandles = [];

phiList = [-pi -4*pi/6:pi/6:0];
pointParams = [];
gammaStops = 3;

for phiIdx = 1:length(phiList)
    thisPhi = phiList(phiIdx);
    gamma = atan2(sin(thisPhi),(1-cos(thisPhi)));
    r1 = r0*(sin(thisPhi)/sin(gamma));
    gammaList = linspace(-gamma,gamma,gammaStops);
    
    if(gamma == 0)
        pointParams = [pointParams; r0,0];
    else
        % convert from r1,gamma back to r2,phi
        for gammaIdx = 1:length(gammaList)
           thisGamma = gammaList(gammaIdx);
           
           phi = atan2(r1*sin(thisGamma),(r0-r1*cos(thisGamma)))
           if(phi == 0 || phi == pi)
                   r2 = abs(r0-r1);
           else
               r2  = r1*sin(thisGamma)/sin(phi);
           end
           pointParams = [pointParams; r2,phi]; 
        end
        
        % plot isovelocity lines
        gammaVals = -1*abs(gamma):0.01:abs(gamma);
        xvals = r1*sin(gammaVals);
        yvals = r1*cos(gammaVals);
        plot(xvals,yvals,'--','color',[0 0.7 0],'LineWidth',2.5);
        
    end
end
%%
for paramIdx = 1:size(pointParams,1)
    
    % get parameters for position of tracked point on wheel
    rho = pointParams(paramIdx,1);
    phi = pointParams(paramIdx,2);
    
    % generate position, velocity, and acceleration
    rx = rho*sin(theta+phi)-r0*theta;
    ry = r0-rho*cos(theta+phi);
    
    % generate only sampled vx and ay
    % ASSUMES CONSTANT ANGULAR VELOCITY AND ANGULAR ACCELERATION
    rx_samp = rho*sin(theta_samp+phi)-r0*theta_samp;
    ry_samp = r0-rho*cos(theta_samp+phi);
    
    % trim trajectory
    trajMaxIdx = find(rx >= trajMaxX,1,'first');
    if( ~isempty(trajMaxIdx))
       rx = rx(1:trajMaxIdx-1);
       ry = ry(1:trajMaxIdx-1);
    end
    
    % show trajectory
    plot(rx,ry,':','color',[0 0 1],'LineWidth',1.5);
    plot(rx(1),ry(1),'.','color',[0 0 1],'MarkerSize',25);
    
    
end

% linkaxes(axHandles,'x');
xlim([-1.1,3]);
axis image
h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';
