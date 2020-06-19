% restart
close all; clear all; clc;

% wheel radius
r0 = 1;

% constant angular velocity
theta_dot = -1; % [rad/s]
theta_ddot = 0; % [rad/s^2]

% time
t=0:pi/16:4*pi;
t_samp = 0:pi/4:4*pi;

% generate angular positions
theta = theta_dot*t;
theta_samp = theta_dot*t_samp;

% plot cycloid path
figure;
set(gcf,'Position',[3.602000e+02 4.180000e+01 8.056000e+02 7.408000e+02]);
hold on; grid on;

vscale = 2;
ascale = 1.5;
axHandles = [];

pointParams = [
    -1.00, 0;
    %-0.75, 0;
    -0.50, 0;
    %-0.25, 0;
    0.00, 0;
    %0.25, 0;
    0.50, 0;
    %0.75, 0;
    1.00, 0];

for paramIdx = 1:size(pointParams,1)
    
    % start a subplot
    axHandles(end+1) = subplot(5,1,paramIdx);
    hold on;
    
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
    
    vx_samp = rho*theta_dot.*cos(theta_samp+phi)-r0*theta_dot;
    vy_samp = rho*theta_dot.*sin(theta_samp+phi);
    
    ax_samp = -rho*theta_dot.^2*sin(theta_samp+phi)+rho*theta_ddot.*cos(theta_samp+phi)-r0*theta_ddot;
    ay_samp = rho*theta_dot.^2*cos(theta_samp+phi)+rho*theta_ddot*sin(theta_samp+phi);
    
    
    % show trajectory
    plot(rx,ry,':','color',[0 0 1],'LineWidth',1.6);
    
    % show velocity and acceleration vectors
    for sampIdx = 1:length(theta_samp)
        thisVX = vx_samp(sampIdx);
        thisVY = vy_samp(sampIdx);
        plot(rx_samp(sampIdx)+[0 vx_samp(sampIdx)/vscale],ry_samp(sampIdx)+[0 vy_samp(sampIdx)/vscale],'-','color',[0 0.7 0],'LineWidth',1.6);
        plot(rx_samp(sampIdx)+[0 ax_samp(sampIdx)/ascale],ry_samp(sampIdx)+[0 ay_samp(sampIdx)/ascale],'-','color',[1 0 0],'LineWidth',1.6);
        plot(rx_samp(sampIdx),ry_samp(sampIdx),'.','color',[0 0 1],'MarkerSize',10);
    end
    
    
    circTheta = 0:0.01:2*pi;
    circX = r0*cos(circTheta);
    circY = r0+r0*sin(circTheta);
    plot(circX,circY,'b-','LineWidth',1.6);
    plot(rho*sin(phi),r0-rho*cos(phi),'.','MarkerSize',25,'Color',[0 0 1]);
    plot(0,r0,'k.','MarkerSize',15);
    axis image
    h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';

end

linkaxes(axHandles,'xy');