% restart
close all; clear; clc;

% true points
p_start = [  1  1  1;
           -1  1  1;
           -1 -1  1;
            1 -1  1;
            1  1 -1;
           -1  1 -1;
           -1 -1 -1;
            1 -1 -1  ];
            
        
plot3(p_start(:,1),p_start(:,2),p_start(:,3),'.','MarkerSize',40,'Color',[0.0 0.0 0.8]);
hold on; grid on; axis equal;
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');

% generate perturbation
% 1. rotate by 15deg about z axis
% 2. transate by [0.1 0.2 0.3]'
theta = 15*pi/180;
t = [0.1 0.2 0.3]';
V = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
T_perturb_t = [eye(3), zeros(3,1); t' 1];
T_perturb_r = [V, zeros(3,1); zeros(1,3) 1];
T_perturb = T_perturb_r*T_perturb_t;   % rotates first, then translates, when premultiplied by a ROW VECTOR
% T_perturb = T_perturb_t*T_perturb_r;   % translates first, then rotates, when premultiplied by a ROW VECTOR
T_perturb

% apply and sho perturbation
p_end = hTF(p_start,T_perturb,1);
plot3(p_end(:,1),p_end(:,2),p_end(:,3),'.','MarkerSize',40,'Color',[0.0 0.8 0.0]);

% we're trying to learn p_perturb from ICP
pc_start = pointCloud(p_start);
pc_end = pointCloud(p_end);
[Tform_icp,movingreg,rmse] = pcregistericp(pc_start,pc_end);  % args: moving, fixed
T_icp = Tform_icp.T % generic, scaled point set frame to "standardized" observed frame

% check error in transformation matrix
T_perturb - T_icp

% show ICP points
p_icp = hTF(p_start,T_icp,1);
plot3(p_icp(:,1),p_icp(:,2),p_icp(:,3),'o','MarkerSize',14,'Color',[0.8 0.0 0.8]);

% finsih plot
legend('Initial Cloud','True Transformation','ICP Transformation'); 
