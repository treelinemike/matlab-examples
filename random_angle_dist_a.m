% restart
close all; clear; clc;

N_samp = 10000;
r_max = 6;

% choose a random angle
theta = 2*pi*rand(N_samp,1);
r = r_max*rand(N_samp,1);
x_samp = [r.*cos(theta) r.*sin(theta)];
xq_vec = -1.2*r_max:0.1:1.2*r_max;
[XX, YY] = meshgrid(xq_vec,xq_vec);
xq = [XX(:) YY(:)];

fq = mvksdensity(x_samp, xq);%,'kernel','epanechnikov','bandwidth',.2);
FQ = reshape(fq,size(XX));

figure;
hold on; grid on; axis equal;
plot(x_samp(:,1),x_samp(:,2),'.','MarkerSize',4,'Color',[0 0 0.8]);
contour(XX,YY,reshape(fq,size(XX)));

figure;
surf(XX,YY,FQ)

figure;
hold on; grid on;
hist(r);
title('\bfDistribution of Sampled Angles');
xlabel('\bfAngle [deg]');
ylabel('\bfCount');