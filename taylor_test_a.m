% restart
close all; clear all; clc;

a = pi;

x = linspace(0,2*pi,1000);
y = sin(x);

y0 = sin(a)*ones(size(x));
y1 = y0 + cos(a)*(x-a);
y2 = y1 - 0.5*sin(a)*(x-a).^2;
y3 = y2 - (1/factorial(3))*cos(a)*(x-a).^3;

figure;
set(gcf,'Position',[0285 0229 8.688000e+02 0420]);
hold on; grid on;
plot(x,y,'b-','LineWidth',3);
plot(x,y0,'r-','LineWidth',1);
plot(x,y1,'-','Color',(1/255)*[244 152 66],'LineWidth',1);
plot(x,y2,'-','Color',(1/255)*[2 2 25],'LineWidth',1);
plot(x,y3,'-','Color',(1/255)*[80 244 66],'LineWidth',1);