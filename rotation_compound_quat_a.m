% test compound rotation with quaternions
% for post-multiplying row vectors

% restart
close all; clear; clc;

% define triad
points = [0 0 0; 1 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 1];

% configure plot and show initial triad position
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

% assemble total rotation
q = quatmult(q2,q1);
R = quat2matrix(q);

% check doing combination in tangent space
ss1 = theta1*[0 -u1(3) u1(2); u1(3) 0 -u1(1);-u1(2) u1(1) 0]';  % note transpose b/c row vectors (ugh)
ss2 = theta2*[0 -u2(3) u2(2); u2(3) 0 -u2(1);-u2(2) u2(1) 0]';  % note transpose b/c row vectors (ugh)
max(max(abs(quat2matrix(q1) - expm(ss1))))
max(max(abs(quat2matrix(q2) - expm(ss2))))
max(max(abs(expm(ss2)*expm(ss1) - R)))

% Test Baker-Campbell-Hausdorff formula for combining matrix exponential
% Not a good idea here because ss1 and ss2 are singular
%max(max(abs( expm(ss2 + ss1 + .5*inv(ss2)*inv(ss1)*ss2*ss1 ) - R)))

% perform total rotation
% note: row vector is post multiplied by rotation matrix
points_new = points*R;

% show rotated triad
plot3(points_new(:,1),points_new(:,2),points_new(:,3),'r-','LineWidth',1.6);