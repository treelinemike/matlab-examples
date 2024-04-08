% Simulate Euler's equations of motion for rotations of an STL file
%
% Requires several physical paramters including the location of the CM and
% the rotation matrix to align principal axes with inertial space
%
% Author:   Mike Kokko
% Modified: 07-Apr-2024

% restart
close all; clear; clc;

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 1;         % [s] simulation end time
dt = 0.005;     % [s] timestep size

% plotting parameters
anim_step = 1;             % skip this many frames to speed up animation
show_translation = false;   % easier to see rotation if translation not shown

% load STL file
% notes:   STL units should be [mm]
%          For Solidworks export, use option "Do not translate STL output data to positive space"
stl_raw = stlread('candycane.stl');

% physical parameters
m   = 1.459;        % [kg] mass
Ixx = 0.001002;     % [kg*m^2] moment of inertia about body x axis
Iyy = 0.004685;     % [kg*m^2] moment of inertia about body y axis
Izz = 0.005523;     % [kg*m^2] moment of inertia about body z axis
Icm  = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];  % inertia matrix taken at CM about principal axes
R_pa = [0.44 0.90 0; -0.90 0.44 0; 0 0 1];
x_cm = [0.024340, 0.108379, 0.0];
params.Icm = Icm;
params.Mx = 0;
params.My = 0;
params.Mz = 0;
params.g = 9.81; % [m/s^2]

% rotate STL file to principal axes
% then view it to make sure the transformation is correct
stl_pts = ((stl_raw.Points/1000)-x_cm)*R_pa';
stl = triangulation(stl_raw.ConnectivityList,stl_pts);
figure;
hold on; grid on; axis equal;
patch('Faces',stl.ConnectivityList,'Vertices',stl.Points,'FaceColor',[0.8 0.2 0.2],'EdgeColor',[0 0 0],'LineWidth',0.5);
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
title('\bfSTL Transformed to Principal Axes about CM');

% initial conditions X0 = [theta_x theta_y theta_z omega_x omega_y omega_z t_x t_y t_z v_x v_y v_z]
X0 = [0 0 0 20 10 0 0 0 0 3*sind(45) 0 3*cosd(45)]'; % [rad rad rad rad/s rad/s rad/s m m m m/s m/s m/s]'
X = X0;

% data storage
time = [t0];
data = [X0];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % propigate state
    [T,X] = ode45(@(t,X) simpleEulerSimStateProp(t,X,params),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

%% Develop and display animation of motion
% define Cartesian frames
triad_scale = 0.7*max(max(abs(stl_pts))); % scaling parameter
triad_XYZ = triad_scale*[0 0 0; 1 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 1];
triad_xyz_template = .8*triad_XYZ;

% keep track of body axes to transform H back
R = eye(3);

% initialize figure
figure;
hold on; grid on;
firstrun = 1;

% animate each frame of results
for tIdx = 2:size(data,2)
    
    % compute rotation between body coordinates and inertial space
    d_theta_x = data(1,tIdx) - data(1,tIdx-1);
    d_theta_y = data(2,tIdx) - data(2,tIdx-1);
    d_theta_z = data(3,tIdx) - data(3,tIdx-1);
    xRot = [1 0 0 ; 0 cos(d_theta_x) -sin(d_theta_x); 0 sin(d_theta_x) cos(d_theta_x)];
    yRot = [cos(d_theta_y) 0 sin(d_theta_y); 0 1 0; -sin(d_theta_y) 0 cos(d_theta_y)];
    zRot = [cos(d_theta_z) -sin(d_theta_z) 0 ; sin(d_theta_z) cos(d_theta_z) 0; 0 0 1];
    R  = R*(xRot*yRot*zRot*eye(3));     % incremental rotations so order shouldn't matter
    if(show_translation)
        t = data(7:9,tIdx);
        t_hist = data(7:9,1:tIdx);
    else
        t = zeros(3,1);
        t_hist = zeros(3,1);
    end

    % skip frames if desired to speed up animation
    % don't do this in the for loop b/c need to update rotation at each step
    if( mod(tIdx-2,anim_step) == 0 )
        % transform point cloud to correct location in inertial space
        triad_xyz = triad_xyz_template*R' + t';
        
        % compute angular velocity in terms of the inertial basis (XYZ)
        omega = data(4:6,tIdx);
        omega_XYZ = R*omega;
        
        % compute angular momentum vector in terms of inertial basis (XYZ)
        % THIS SHOULD STAY CONSTANT (no external moments)
        Hcm_xyz = Icm*omega;
        Hcm_XYZ = R*Hcm_xyz;
        
        % clear axes and start plotting the current frame
        cla;
        
        % plot trajectory
        plot3(t_hist(1,:),t_hist(2,:),t_hist(3,:),':','LineWidth',2,'Color',[0.8 0 0]);

        % plot XYZ (all black) and xyz (x=red, y=green, z=blue) coordinate frames
        plot3(t(1)+triad_XYZ(:,1),t(2)+triad_XYZ(:,2),t(3)+triad_XYZ(:,3),'LineWidth',4,'Color','k')
        plot3([triad_xyz(1,1) triad_xyz(2,1)],[triad_xyz(1,2) triad_xyz(2,2)],[triad_xyz(1,3) triad_xyz(2,3)],'LineWidth',3,'Color',[0.7 0 0]);
        plot3([triad_xyz(3,1) triad_xyz(4,1)],[triad_xyz(3,2) triad_xyz(4,2)],[triad_xyz(3,3) triad_xyz(4,3)],'LineWidth',3,'Color',[0 0.7 0]);
        plot3([triad_xyz(5,1) triad_xyz(6,1)],[triad_xyz(5,2) triad_xyz(6,2)],[triad_xyz(5,3) triad_xyz(6,3)],'LineWidth',3,'Color',[0 0 0.7]);
        
        % normalize and plot angular velocity and momentum
        omega_norm = 2.6*omega_XYZ/norm(omega_XYZ);
        Hcm_norm = 2.6*Hcm_XYZ/norm(Hcm_XYZ);
        ph(1) = plot3(t(1)+[0 omega_norm(1)],t(2)+[0 omega_norm(2)],t(3)+[0 omega_norm(3)],':','LineWidth',3','Color',[1 0 1]);
        ph(2) = plot3(t(1)+[0 Hcm_norm(1)],t(2)+[0 Hcm_norm(2)],t(3)+[0 Hcm_norm(3)],':','LineWidth',3','Color',[0 1 1]);
        
        % plot patch object
        patch('Faces',stl.ConnectivityList,'Vertices',stl_pts*R' + t','FaceColor',[0.8 0.2 0.2],'EdgeColor',[0 0 0],'LineWidth',0.5);
        
        % finish formatting axes
        axis equal;
        if(show_translation)
            xlim([-sqrt(2)*triad_scale 1]);
            ylim([-sqrt(2)*triad_scale sqrt(2)*triad_scale]);
            zlim([-0.1 0.3]);
        else
            xlim(sqrt(2)*[-triad_scale triad_scale]);
            ylim(sqrt(2)*[-triad_scale triad_scale]);
            zlim(sqrt(2)*[-triad_scale triad_scale]);
        end

        xlabel('\bfx');
        ylabel('\bfy');
        zlabel('\bfz');
        
        % add legend on first pass
        % except legends don't work so well on Mac for some reason?
        if(firstrun)
            %         legend(ph,{'Angular Velocity','Angular Momentum'},'Location','southoutside','AutoUpdate','off');
            firstrun = 0;
        end
        view([145,30]);
        
        drawnow;
    end
end

% propagate state
function Xdot = simpleEulerSimStateProp(t,X,params)

% recover parameters
Icm = params.Icm;
Ixx = Icm(1,1);
Iyy = Icm(2,2);
Izz = Icm(3,3);
Mx = params.Mx;
My = params.My;
Mz = params.Mz;

% deconstruct state vector
theta_x    = X(1);
theta_y    = X(2);
theta_z    = X(3);
omega_x    = X(4);
omega_y    = X(5);
omega_z    = X(6);
t_x        = X(7);
t_y        = X(8);
t_z        = X(9);
v_x        = X(10);
v_y        = X(11);
v_z        = X(12);

% construct Xdot from differential equation
% note:     X    = [theta_x theta_y theta_z omega_x omega_y omega_z t_x t_y t_z v_x v_y v_z]
% therefore Xdot = [omega_x omega_y omega_z alpha_x alpha_y alpha_z v_x v_y v_z a_x a_y a_z]
Xdot = zeros(12,1);
Xdot(1,:) = omega_x;
Xdot(2,:) = omega_y;
Xdot(3,:) = omega_z;
Xdot(4,:) = (Mx - omega_y*omega_z*(Izz-Iyy))/Ixx;
Xdot(5,:) = (My - omega_x*omega_z*(Ixx-Izz))/Iyy;
Xdot(6,:) = (Mz - omega_x*omega_y*(Iyy-Ixx))/Izz;
Xdot(7,:) = v_x;
Xdot(8,:) = v_y;
Xdot(9,:) = v_z;
Xdot(10,:) = 0;
Xdot(11,:) = 0;
Xdot(12,:) = -1*params.g;


end