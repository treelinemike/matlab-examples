% restart
close all; clear all; clc;

theta = 0:0.001:pi;
err = (theta-sin(theta))./theta *100;

figure;
set(gcf,'Position',[1.426000e+02 5.058000e+02 1.348800e+03 2.560000e+02]);
hold on; grid on;
plot(theta*180/pi,err,'-','LineWidth',1.6);
xlabel('\bfAngle [deg]');
ylabel('\bfLinerization Error [%]');
