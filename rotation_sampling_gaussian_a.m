
% Test sampling rotations from a Gaussian
% following: http://ethaneade.com/lie_groups.pdf
%
% Essentially sample from a multivariate gaussian in so(3) which is the
% tangent space of SO(3) about the identity eye(3)...
% These 3-vectors v will ALWAYS* be valid rotations when transformed via
% exponential map: R = expm([0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]);
% * may need to prove this?

% restart
close all; clear; clc;
rng('default');

% options
N_samp = 1000;
v_nom = [1 0 0];  % nominal vector that we are going to perturb
% ts_mean = [0 -pi/2 0]';
ts_mean = [ 0   -1.5459   1.5459]'; % tangent space coefficients of mean perturbation rotation
                                    % note: NOT how much to rotate about
                                    % each axis; rather: angle*[u1 u2 u3]'
sigma_x = (5)*pi/180;  % affects x component of rotation axis unit vector
sigma_y = (5)*pi/180;  % affects y component of rotation axis unit vector
sigma_z = (5)*pi/180;  % affects z component of rotation axis unit vector
k_tang = 0.3;

% data storage
v_rot = zeros(3,N_samp);

% define covariance in tangent space
% that is, vectors whose elements are [theta_x, theta_y, theta_z]
COV_tang = diag([(sigma_x)^2, (sigma_y)^2, (sigma_z)^2]);
samp = mvnrnd(ts_mean',COV_tang,N_samp);  % gaussian sampling
samp_angax = zeros(4,N_samp);
% samp = 2*(pi/2)*rand(N_samp,3)-(pi/2); % test uniform sampling instead

% not sure if this is correct: deterministic uniform mapping
% % unifSpace = -pi/2:0.1:pi/2;
% % [X,Y,Z] = meshgrid(unifSpace,unifSpace,unifSpace);
% % samp = [X(:) Y(:) Z(:)];

% prepare plot
figure;
set(gcf,'Position',[0378 0393 1225 0460]);
subplot(1,2,1);
hold on; grid on; axis equal;
title('\bfSO(3) Samples Applied to Vector in R^3');
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
view([45 36]);

% convert samples to rotations
for sampIdx = 1:size(samp,1)

   % convert sample from so(2) coefficients to quaternion 
   q_samp = tang2quat(samp(sampIdx,:)');
   
   % apply quaternion rotation to nominal vector
   v_rot(:,sampIdx) = quatrotate(q_samp,v_nom);
   
   % convert to angle/axis and plot AXIS OF ROTATION for this sample
   angax = quat2angax(q_samp);
   tang = quat2tang(q_samp);
   samp_angax(:,sampIdx) = angax;
%    plot3([0 angax(2)],[0 angax(3)],[0 angax(4)],'-','LineWidth',0.1,'Color',[0.8 0 0.8]);
   plot3(k_tang*[0 tang(1)],k_tang*[0 tang(2)],k_tang*[0 tang(3)],'-','LineWidth',0.1,'Color',[0.8 0 0.8]);
   
end

% show random rotations of the nominal vector
plot3(v_rot(1,:),v_rot(2,:),v_rot(3,:),'.','MarkerSize',2,'Color',[0 0 0.8]);

% draw references
theta = 0:0.01:2*pi;
circxy = [cos(theta); sin(theta); zeros(size(theta))];
circxz = [cos(theta); zeros(size(theta)); sin(theta)];
circyz = [zeros(size(theta)); cos(theta); sin(theta)];
plot3(circxy(1,:),circxy(2,:),circxy(3,:),'-','LineWidth',1.6,'Color',[0.8 0 0]);
plot3(circxz(1,:),circxz(2,:),circxz(3,:),'-','LineWidth',1.6,'Color',[0.8 0 0]);
plot3(circyz(1,:),circyz(2,:),circyz(3,:),'-','LineWidth',1.6,'Color',[0.8 0 0]);
plot3(0,0,0,'.','MarkerSize',25,'Color',[0 0 0.8]);

% original (pre-rotated) vector that random rotations are operating on
plot3([0 v_nom(1)],[0 v_nom(2)],[0 v_nom(3)],'--','LineWidth',1.6,'Color',[0 0 0.8]);

% find mean of samples in tangent space
% then compute and apply corresponding rotation
samp_mean = mean(samp,1)';
samp_mean_vec = quatrotate(tang2quat(samp_mean),v_nom);
plot3([0 samp_mean_vec(1)],[0 samp_mean_vec(2)],[0 samp_mean_vec(3)],'-','LineWidth',1.6,'Color',[0 0 0.8]);
plot3(samp_mean_vec(1),samp_mean_vec(2),samp_mean_vec(3),'.','MarkerSize',20,'Color',[0 0 0]);

% show the desired mean axis of rotation specified in the gaussian
% distribution that we sample from
ts_mean_uv = unitvec(ts_mean);
samp_mean_uv = unitvec(samp_mean);
plot3([0 ts_mean_uv(1)],[0 ts_mean_uv(2)],[0 ts_mean_uv(3)],'-','Color',[0 0.8 0]);
plot3(k_tang*ts_mean(1),k_tang*ts_mean(2),k_tang*ts_mean(3),'.','MarkerSize',20,'Color',[0 0.8 0]);
plot3([0 samp_mean_uv(1)],[0 samp_mean_uv(2)],[0 samp_mean_uv(3)],'-','Color',[0 0 0]);
plot3(k_tang*samp_mean(1),k_tang*samp_mean(2),k_tang*samp_mean(3),'.','MarkerSize',20,'Color',[0 0 0]);

% show tangent space
subplot(1,2,2);
hold on; grid on; axis equal;
title('\bfCoefficients of Generators in Tangent Space so(3)');
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
view([45 36]);
plot3(samp(:,1),samp(:,2),samp(:,3),'.','MarkerSize',2,'Color',[0.8 0 0.8]);
plot3(ts_mean(1),ts_mean(2),ts_mean(3),'.','MarkerSize',20,'Color',[0 0.8 0]);
plot3(samp_mean(1),samp_mean(2),samp_mean(3),'.','MarkerSize',20,'Color',[0 0 0]);
legend('Samples','Desired Mean','Sample Mean','Location','NorthEast');

% display results
fprintf('    Desired >> Theta: %+8.4f rad; Axis: [%+8.4f,%+8.4f,%+8.4f ]\n',tang2angax(ts_mean));
fprintf('Sample Mean >> Theta: %+8.4f rad; Axis: [%+8.4f,%+8.4f,%+8.4f ]\n',tang2angax(samp_mean));
