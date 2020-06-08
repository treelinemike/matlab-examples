% Test sampling rotations from a Uniform or Gaussian distribution
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

% general options
% SAMPLE_TYPE = 'uniform';   % 'uniform', 'uniform_angle', or 'gaussian'
SAMPLE_TYPE = 'uniform_angle';
% SAMPLE_TYPE = 'gaussian';
N_samp = 10000;
v_nom = [1 0 0]';  % nominal vector that we are going to perturb
k_tang = 0.3; % scale length of rotation axes for display purposes only
% ts_mean = [0 0 0]';
ts_mean = [ 0   -1.5459   1.5459]'; % tangent space coefficients of mean perturbation rotation
                                    % note: NOT how much to rotate about
                                    % each axis; rather: angle*[u1 u2 u3]'

                                    
                                    
% UNIFORM OPTIONS
theta_max = 15*pi/180;

% GAUSSIAN OPTIONS
sigma_x = (5)*pi/180;  % affects x component of rotation axis unit vector
sigma_y = (5)*pi/180;  % affects y component of rotation axis unit vector
sigma_z = (5)*pi/180;  % affects z component of rotation axis unit vector

% data storage
samp_angax = zeros(4,N_samp);

% assemble std. dev. vector
sigma = [sigma_x, sigma_y, sigma_z];

% samples in tangent space
switch SAMPLE_TYPE
    case 'uniform'
        [samp,delta_theta] = randRotUnif_t(ts_mean,theta_max,N_samp);
    case 'uniform_angle'
        [samp,delta_theta] = randRotUnifAng_t(ts_mean,theta_max,N_samp);    
    case 'gaussian'
        [samp,delta_theta] = randRotGauss_t(ts_mean,sigma,N_samp);
    otherwise
        error('Unknown sampling type. Choose ''uniform'' or ''gaussian''');
end
% rotated vectors
v_rot = rot1vec_t(v_nom,samp);

% prepare plot
figure;
set(gcf,'Position',[1.626000e+02 2.266000e+02 1.225600e+03 0468]);
subplot(1,2,1);
hold on; grid on; axis equal;
title('\bfSO(3) Samples Applied to Vector in R^3');
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
view([45 36]);

% display axes of rotation, with magnitudes scaled to represent rotation
% angles
for sampIdx = 1:N_samp

    % extract net rotation and convert to quaternion
    t_total = samp(:,sampIdx);
    q_total = tang2quat(t_total);
    
    % convert to angle/axis and plot AXIS OF ROTATION for this sample
    angax = quat2angax(q_total);
    samp_angax(:,sampIdx) = angax;
    
    %    plot3([0 angax(2)],[0 angax(3)],[0 angax(4)],'-','LineWidth',0.1,'Color',[0.8 0 0.8]);
    plot3(k_tang*[0 t_total(1)],k_tang*[0 t_total(2)],k_tang*[0 t_total(3)],'-','LineWidth',0.1,'Color',[0.8 0 0.8]);
   
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
samp_mean = mean(samp,2);
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
plot3(samp(1,:),samp(2,:),samp(3,:),'.','MarkerSize',4,'Color',[0.8 0 0.8]);
plot3(ts_mean(1),ts_mean(2),ts_mean(3),'.','MarkerSize',20,'Color',[0 0.8 0]);
plot3(samp_mean(1),samp_mean(2),samp_mean(3),'.','MarkerSize',20,'Color',[0 0 0]);
legend('Samples','Desired Mean','Sample Mean','Location','NorthEast');

% display results
fprintf('    Desired >> Theta: %+8.4f rad; Axis: [%+8.4f,%+8.4f,%+8.4f ]\n',tang2angax(ts_mean));
fprintf('Sample Mean >> Theta: %+8.4f rad; Axis: [%+8.4f,%+8.4f,%+8.4f ]\n',tang2angax(samp_mean));
fprintf('Max perturbation angle: %6.2f rad = %6.2f°\n',max(delta_theta)*[1 180/pi]);

% Show distribution of angles
% NOTE: This will NOT be uniform for a uniform distribution
% because we sample in tangent space, and the probability of samping a
% region is proportional to the volume of that region. In the sampled
% sphere space there is more volume associated with larger angles. 
% See excel file for theoretical comparison, doesn't quite line up
% (histogram goes with ang^2, theory with ang^3.
% Consider interplay between inclination of axis and angle of rotation?
figure;
hold on; grid on;
hist(delta_theta*180/pi);
title('\bfDistribution of Sampled Angles');
xlabel('\bfAngle [deg]');
ylabel('\bfCount');
