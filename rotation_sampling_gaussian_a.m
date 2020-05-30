% Test sampling rotations from a gaussian
% following: http://ethaneade.com/lie_groups.pdf

% restart
close all; clear; clc;

% options
sigma_x = (10)*pi/180;
sigma_y = (5)*pi/180;
sigma_z = (1)*pi/180;
N_samp = 10000;
v_nom = [0 0 1];  % nominal vector that we are going to perturb

% data storage
v_samp = zeros(3,N_samp);

% define covariance in tangent space
% that is, vectors whose elements are [theta_x, theta_y, theta_z]
COV_tang = diag([(sigma_x)^2, (sigma_y)^2, (sigma_z)^2]);
samp = mvnrnd(zeros(1,3),COV_tang,N_samp);

% convert samples to rotations
for sampIdx = 1:size(samp,1)
   ts_vec = samp(sampIdx,:)';  % vector in tangent space (coefficients of generators)
   ts_ssmat = [0 -ts_vec(3) ts_vec(2); ts_vec(3) 0 -ts_vec(1); -ts_vec(2) ts_vec(1) 0];  % skew-symmetric (cross product) matrix from tangent space
   R_samp = expm(ts_ssmat);  % exponential map transforms element of algebra so(3) into element of group SO(3)
   assert( abs(det(R_samp)-1) < 1e-6,'Invalid rotation matrix!')
   q_samp = matrix2quat(R_samp);

   v_samp(:,sampIdx) = quatrotate(q_samp,v_nom);
   
end

% prepare plot
figure;
hold on; grid on; axis equal;
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
view([45 36]);
plot3(v_samp(1,:),v_samp(2,:),v_samp(3,:),'.','MarkerSize',5,'Color',[0.8 0 0.8]);

% draw references
theta = 0:0.01:2*pi;
circxy = [cos(theta); sin(theta); zeros(size(theta))];
circxz = [cos(theta); zeros(size(theta)); sin(theta)];
circyz = [zeros(size(theta)); cos(theta); sin(theta)];
plot3(circxy(1,:),circxy(2,:),circxy(3,:),'-','LineWidth',1.6,'Color',[0.8 0 0]);
plot3(circxz(1,:),circxz(2,:),circxz(3,:),'-','LineWidth',1.6,'Color',[0.8 0 0]);
plot3(circyz(1,:),circyz(2,:),circyz(3,:),'-','LineWidth',1.6,'Color',[0.8 0 0]);
plot3(0,0,0,'.','MarkerSize',25,'Color',[0 0 0.8]);
plot3([0 v_nom(1)],[0 v_nom(2)],[0 v_nom(3)],'-','LineWidth',1.6,'Color',[0 0 0.8]);
title('\bfGaussian SO(3) Samples Applied to Vector in R^3');
drawnow;