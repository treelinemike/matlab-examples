% restart
close all; clear; clc;

% original vectors
a_colvec = [3 4]';
a_rowvec = a_colvec';

% rotation
theta = 20*pi/180;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% rotated vectors
a_prime_colvec = R*a_colvec;
a_prime_rowvec = a_rowvec*R';

% plot results
figure;
set(gcf,'Position',[0417 3.882000e+02 0560 3.152000e+02]);
subplot(1,2,1);
hold on; grid on; axis equal;
xlim([-10 10]);
ylim([-10 10]);
plot([0 a_colvec(1)],[0 a_colvec(2)],'b-','LineWidth',1.);
plot([0 a_prime_colvec(1)],[0 a_prime_colvec(2)],'r-','LineWidth',1.);
legend('Original Vector','Rotated Vector');
title('\bfColumn Vectors');

subplot(1,2,2);
hold on; grid on; axis equal;
xlim([-10 10]);
ylim([-10 10]);
plot([0 a_rowvec(1)],[0 a_rowvec(2)],'b-','LineWidth',1.);
plot([0 a_prime_rowvec(1)],[0 a_prime_rowvec(2)],'r-','LineWidth',1.);
legend('Original Vector','Rotated Vector');
title('\bfRow Vectors');