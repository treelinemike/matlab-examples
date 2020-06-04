
% Test sampling rotations from a Uniform distribution on the 3-ball
% (interior of S^2)
% following: http://ethaneade.com/lie_groups.pdf
% and also: http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
%
% Essentially take three samples from U(-1 1); reject if not in the 3-ball, and multiply by the desired max rotation angle
% The results is an element of the tangent space of SO(3) [ so(3) ] about the identity eye(3)...
% These 3-vectors v will ALWAYS* be valid rotations when transformed via
% exponential map: R = expm([0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]);
% * may need to prove this?

% restart
close all; clear; clc;
rng('default');

% options
N_samp = 4000;
v_nom = [1 0 0];  % nominal vector that we are going to perturb
% ts_mean = [ 0   0   0 ]';
ts_mean = [ 0   -1.5459   1.5459]'; % tangent space coefficients of mean perturbation rotation
                                    % note: NOT how much to rotate about
                                    % each axis; rather: angle*[u1 u2 u3]'
theta_max = 15*pi/180;
k_tang = 0.3;

% get mean rotation
q_mean = tang2quat(ts_mean);

% data storage
v_rot = zeros(3,N_samp);
samp = zeros(N_samp,3);
samp_angax = zeros(4,N_samp);

% samplue using rejection method
N_accepted = 0;
loopCount = 0;

while( N_accepted < N_samp )
    U_samp = (2*rand(1,3)-1);
    % Note: Typically we would need to draw another random sample here,
    % this time from u(0,1) and compare that value to the ratio of the
    % f(x)/(c*g(x)) where f(x) is the density we want to sample from and
    % c*g(x) is the envelope function. Here we don't need the random draw
    % because f(x)/(c*g(x)) = 1 for all x inside the unit sphere and 0
    % outside. Thus, samples are always accepted if their norm is less than
    % 1.0.
    if( norm(U_samp) <= 1 )
        N_accepted = N_accepted + 1;
        samp(N_accepted,:) =  quat2tang(quatmult(q_mean, tang2quat(theta_max*U_samp) )) ;
    end
    loopCount = loopCount + 1;
end
fprintf('Acceptance ratio: %5.2f%%\n',100*N_samp/loopCount);





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
plot3(samp(:,1),samp(:,2),samp(:,3),'.','MarkerSize',4,'Color',[0.8 0 0.8]);
plot3(ts_mean(1),ts_mean(2),ts_mean(3),'.','MarkerSize',20,'Color',[0 0.8 0]);
plot3(samp_mean(1),samp_mean(2),samp_mean(3),'.','MarkerSize',20,'Color',[0 0 0]);
legend('Samples','Desired Mean','Sample Mean','Location','NorthEast');

% display results
fprintf('    Desired >> Theta: %+8.4f rad; Axis: [%+8.4f,%+8.4f,%+8.4f ]\n',tang2angax(ts_mean));
fprintf('Sample Mean >> Theta: %+8.4f rad; Axis: [%+8.4f,%+8.4f,%+8.4f ]\n',tang2angax(samp_mean));
