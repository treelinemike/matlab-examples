% Test duality of rotations in tangent space: i.e. so(3)
%
% Tangent space represntation is not unique: for any rotation you can
% equivalently rotate by 2*pi-theta about the same axis, and the resulting
% tangent space representation is on the same line as the original through
% the origin, but not as easy to distinguish as with quaternions where the
% two representations simply have opposite signs. Here we check that for
% small(ish) rotations - those whose angle is < 180deg) we can always
% choose the rotation with smaller angle of rotation (<180deg) to maintain
% a unique representation that can be sampled, averaged, etc...
%
% Max Rotation Angle          Dual Rotation Angle           Margin
% ------------------          --------------------          ------
%       10 deg                       350 deg                340 deg
%       30 deg                       330 deg                300 deg
%       60 deg                       300 deg                240 deg
%       90 deg                       270 deg                180 deg
%      120 deg                       240 deg                120 deg
%      180 deg                       180 deg                  0 deg

% restart
close all; clear all; clc;
rng('default');

% options
N_samp = 2000;
theta_max = 60*pi/180;

% data storage, etc.
loop_counter = 0;
N_accepted = 0;
t_samp = zeros(3,N_samp);

% sample rotations uniformly via rejection method
while N_accepted < N_samp  
    x = (2*rand(3,1)-1);
    if( norm(x) <= 1.0 )
        N_accepted = N_accepted + 1;
        t_samp(:,N_accepted) = theta_max*x;
    end
    loop_counter = loop_counter + 1;
end
fprintf('Rejction method sampling complete, acceptance ratio: %6.2f%%\n',100*N_accepted/loop_counter);

% show samples
figure;
hold on; grid on; axis equal;
plot3(t_samp(1,:),t_samp(2,:),t_samp(3,:),'.','MarkerSize',4,'Color',[0 0 0]);

% prepare for  
neg_count = 0;
maxmin = [0 inf];
for sampIdx = 1:N_samp
   
    t = t_samp(:,sampIdx);
    theta = zeros(2,1);
    theta(1) = norm(t);          % this one will be positive
    theta(2) = theta(1) - 2*pi;  % this one will be negative
    
    % sort angles so the angle closeset to zero is first
    % this should only be necessary if theta_max > 180deg
    [~,so] = sort(abs(theta));
    if(so(1) ~= 1)
        warning('Corrected sort order.');
    end
    theta = theta(so);
    
    % extract unit vector for this rotation
    u1 = t/theta(1);
    
    % compute equivalent tangent space vectors
    t_equiv = [ theta(1)*u1 theta(2)*u1];
    q_equiv = tang2quat(t_equiv);
    if(sign(q_equiv(1,1)) == -1)
        neg_count = neg_count + 1
    end
    
    % compute norms: distance of each point from origin
    norms = vecnorm(t_equiv);
    if( norms(1) > maxmin(1) )
        maxmin(1) = norms(1);
    end
    if( norms(2) < maxmin(2) )
        maxmin(2) = norms(2);
    end
    
    % plot results
    plot3(t_equiv(1,:),t_equiv(2,:),t_equiv(3,:),'-','LineWidth',1,'Color',[0 0 0]);
    plot3(t_equiv(1,1),t_equiv(2,1),t_equiv(3,1),'.','MarkerSize',20,'Color',[0.8 0 0]);
    plot3(t_equiv(1,2),t_equiv(2,2),t_equiv(3,2),'.','MarkerSize',20,'Color',[0 0 0.8]);
end
fprintf('Radial disp. in so(3) (i.e. rotation margin) btwn primary and dual reps: %6.2f deg\n',(maxmin(2)-maxmin(1))*180/pi);