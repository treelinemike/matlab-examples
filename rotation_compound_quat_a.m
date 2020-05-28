% restart
close all; clear; clc;

points = [0 0 0; 1 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 1];

figure;
hold on; grid on; axis equal;
plot3(points(:,1),points(:,2),points(:,3),'k-','LineWidth',1.6);
xlim([-3 3]);
ylim([-3 3]);
zlim([-3 3]);
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
view([60 30]);

% primary rotation: +90deg about z axis
theta1 = 90*pi/180;
u1 = [0 0 1]';
q1 = [cos(theta1/2); -sin(theta1/2)*u1]; % note negative sign b/c we'll be applying matrix to ROW VECTORS

% perturbation rotation: +10deg about y axis
theta2 = 10*pi/180;
u2 = [0 1 0]';
q2 = [cos(theta2/2); -sin(theta2/2)*u2]; % note negative sign b/c we'll be applying matrix to ROW VECTORS

q = quatmult(q2,q1);

points_new = points*quat2matrix(q);

plot3(points_new(:,1),points_new(:,2),points_new(:,3),'r-','LineWidth',1.6);