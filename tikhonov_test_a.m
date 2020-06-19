% restart
close all; clear all; clc;

sigma = 1;

x_truth = [0.5; 3.7];
a_truth = (0:0.1:10)';
b_truth = x_truth(1) + x_truth(2)*a_truth + sigma*randn(size(a_truth));
data = [a_truth b_truth];

A = [ones(size(data,1),1) data(:,1)];
b = data(:,2);

x_lsq = (A'*A)\(A'*b)

G = eye(2);  %%% CHOICE HERE! %%%
x_tik = (A'*A + G'*G)\(A'*b)

figure;
hold on; grid on;
plot(data(:,1),data(:,2),'b.','MarkerSize',25);
x = data(:,1);
y_lsq = x_lsq(1) + x_lsq(2)*data(:,1);
y_tik = x_tik(1) + x_tik(2)*data(:,1);

plot(x,y_lsq,'-','Color',[0 0.7 0]);
plot(x,y_tik,'m-');
