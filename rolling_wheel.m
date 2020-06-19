% restart
close all; clear all; clc;

r0 = 1;

% generate trajectory
t=0:pi/4:4*pi;
theta_dot = -1; % [rad/s]
theta = theta_dot*t;

rx = r0*sin(theta)-r0*theta;
ry = r0-r0*cos(theta);

vx = r0*theta_dot*cos(theta)-r0*theta_dot;
vy = r0*theta_dot*sin(theta);

ax = -r0*(theta_dot.^2).*sin(theta);
ay = r0*(theta_dot.^2).*cos(theta);

% plot cycloid path
figure;
set(gcf,'Position',[3.602000e+02 2.674000e+02 8.096000e+02 3.968000e+02]);
plot(rx,ry,'b-','LineWidth',1.6);
hold on; grid on;

vscale = 2;
ascale = 1.5;

for vidx = 1:length(vx)
   x0 = rx(vidx);
   y0 = ry(vidx);
   xvf = x0+vx(vidx)/vscale;
   yvf = y0+vy(vidx)/vscale;
   xaf = x0+ax(vidx)/ascale;
   yaf = y0+ay(vidx)/ascale;
   
   plot([x0 xvf],[y0 yvf],'color',[0 0.7 0],'LineWidth',1.6);
   plot([x0 xaf],[y0 yaf],'r','LineWidth',1.6);
   
end

axis image
title('\bfRolling Wheel Trajectory, Velocity, and Accelration','FontSize',12);
set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});