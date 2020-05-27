% restart
close all; clear; clc;

% original vectors
a_colvec = [3 4 0]';
a_rowvec = a_colvec';

% rotation
theta = 20*pi/180;
R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0 ; 0 0 1];
u = [ 0 0 1 ]';
q = [cos(theta/2); sin(theta/2)*u];

% rotated vectors
a_prime_colvec = R*a_colvec;
a_prime_rowvec = a_rowvec*R';
a_prime_colvec_q = quatrotate(q,a_colvec);
a_prime_rowvec_q = 1; % TODO: FIX...

% plot results
figure;
set(gcf,'Position',[0417 0097 8.712000e+02 6.064000e+02]);

subplot(2,2,1);
hold on; grid on; axis equal;
xlim([-10 10]);
ylim([-10 10]);
plot3([0 a_colvec(1)],[0 a_colvec(2)],[0 a_colvec(3)],'b-','LineWidth',1.);
plot3([0 a_prime_colvec(1)],[0 a_prime_colvec(2)],[0 a_prime_colvec(3)],'r-','LineWidth',1.);
legend('Original Vector','Rotated Vector');
title('\bfColumn Vectors - Matrix');

subplot(2,2,2);
hold on; grid on; axis equal;
xlim([-10 10]);
ylim([-10 10]);
plot3([0 a_rowvec(1)],[0 a_rowvec(2)],[0 a_rowvec(3)],'b-','LineWidth',1.);
plot3([0 a_prime_rowvec(1)],[0 a_prime_rowvec(2)],[0 a_prime_rowvec(3)],'r-','LineWidth',1.);
legend('Original Vector','Rotated Vector');
title('\bfRow Vectors - Matrix');

subplot(2,2,3);
hold on; grid on; axis equal;
xlim([-10 10]);
ylim([-10 10]);
plot3([0 a_colvec(1)],[0 a_colvec(2)],[0 a_colvec(3)],'b-','LineWidth',1.);
plot3([0 a_prime_colvec_q(1)],[0 a_prime_colvec_q(2)],[0 a_prime_colvec_q(3)],'r-','LineWidth',1.);
legend('Original Vector','Rotated Vector');
title('\bfColumn Vectors - Quaternion');

subplot(2,2,4);
hold on; grid on; axis equal;
xlim([-10 10]);
ylim([-10 10]);
plot3([0 a_rowvec(1)],[0 a_rowvec(2)],[0 a_rowvec(3)],'b-','LineWidth',1.);
plot3([0 a_prime_rowvec(1)],[0 a_prime_rowvec(2)],[0 a_prime_rowvec(3)],'r-','LineWidth',1.);
legend('Original Vector','Rotated Vector');
title('\bfRow Vectors - Quaternion');