% Simulate Euler's equations of motion for rotations of an STL file
%
% Requires several physical paramters including the location of the CM and
% the rotation matrix to align principal axes with inertial space
%
% Author:   Mike Kokko
% Modified: 17-Dec-2020

% restart
close all; clear; clc;

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 15;        % [s] simulation end time
dt = 0.005;     % [s] timestep size

% plotting parameters
anim_step = 10;      % skip this many frames to speed up animation
z_shadow = -0.1;     % z coordinate of plane for plotting shadow
h_resample = 0.001;  % discritization for sampling shadow, smaller = slower & less precise

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

% initial conditions X0 = [theta_x_0 theta_y_0 theta_z_0 omega_x_0 omega_y_0 omega_z_0]
X0 = [0 0 0 0 2 0]'; % [rad rad rad rad/s rad/s rad/s]'
X = X0;

% data storage
time = [t0];
data = [X0];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % propigate state
    [T,X] = ode45(@(t,X) simpleEulerSimStateProp(t,X,Icm),odeTime,X);
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
    
    % skip frames if desired to speed up animation
    % don't do this in the for loop b/c need to update rotation at each step
    if( mod(tIdx-2,anim_step) == 0 )
        % transform point cloud to correct location in inertial space
        triad_xyz = triad_xyz_template*R';
        
        % compute angular velocity in terms of the inertial basis (XYZ)
        omega = data(4:6,tIdx);
        omega_XYZ = R*omega;
        
        % compute angular momentum vector in terms of inertial basis (XYZ)
        % THIS SHOULD STAY CONSTANT (no external moments)
        Hcm_xyz = Icm*omega;
        Hcm_XYZ = R*Hcm_xyz;
        
        % clear axes and start plotting the current frame
        cla;
        
        % plot XYZ (all black) and xyz (x=red, y=green, z=blue) coordinate frames
        plot3(triad_XYZ(:,1),triad_XYZ(:,2),triad_XYZ(:,3),'LineWidth',4,'Color','k')
        plot3([triad_xyz(1,1) triad_xyz(2,1)],[triad_xyz(1,2) triad_xyz(2,2)],[triad_xyz(1,3) triad_xyz(2,3)],'LineWidth',3,'Color',[0.7 0 0]);
        plot3([triad_xyz(3,1) triad_xyz(4,1)],[triad_xyz(3,2) triad_xyz(4,2)],[triad_xyz(3,3) triad_xyz(4,3)],'LineWidth',3,'Color',[0 0.7 0]);
        plot3([triad_xyz(5,1) triad_xyz(6,1)],[triad_xyz(5,2) triad_xyz(6,2)],[triad_xyz(5,3) triad_xyz(6,3)],'LineWidth',3,'Color',[0 0 0.7]);
        
        % normalize and plot angular velocity and momentum
        omega_norm = 2.6*omega_XYZ/norm(omega_XYZ);
        Hcm_norm = 2.6*Hcm_XYZ/norm(Hcm_XYZ);
        ph(1) = plot3([0 omega_norm(1)],[0 omega_norm(2)],[0 omega_norm(3)],':','LineWidth',3','Color',[1 0 1]);
        ph(2) = plot3([0 Hcm_norm(1)],[0 Hcm_norm(2)],[0 Hcm_norm(3)],':','LineWidth',3','Color',[0 1 1]);
        
        % plot patch object
        patch('Faces',stl.ConnectivityList,'Vertices',stl_pts*R','FaceColor',[0.8 0.2 0.2],'EdgeColor',[0 0 0],'LineWidth',0.5);
        
        % compute and plot shadow (with its properties)
        newtri = triangulation(stl.ConnectivityList,stl_pts*R');
        poly_shadow = tri_z_proj(newtri,h_resample);
        [xc, yc] = centroid(poly_shadow);
        A = area(poly_shadow);
        fprintf('Polygon area: %6.2e\n',A);
        fprintf('Polygon centroid: (%6.2e,%6.2e)\n',xc,yc);
        t = hgtransform('Matrix',[ eye(3) [0 0 z_shadow]'; [0 0 0 1]]); 
        ph = plot(poly_shadow,'Parent',t);
        plot3(0,0,z_shadow,'.','MarkerSize',20,'Color',[0 0 0]);
        plot3(xc,yc,z_shadow,'.','MarkerSize',20,'Color',[0.8 0 0]);
        
        % finish formatting axes
        axis equal;
        xlim(sqrt(2)*[-triad_scale triad_scale]);
        ylim(sqrt(2)*[-triad_scale triad_scale]);
        zlim(sqrt(2)*[-triad_scale triad_scale]);
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
function Xdot = simpleEulerSimStateProp(t,X,Icm)

% recover moments of inertia
Ixx = Icm(1,1);
Iyy = Icm(2,2);
Izz = Icm(3,3);

% deconstruct state vector
theta_x    = X(1);
theta_y    = X(2);
theta_z    = X(3);
omega_x    = X(4);
omega_y    = X(5);
omega_z    = X(6);

% construct Xdot from differential equation
% note:     X    = [theta_x theta_y theta_z omega_x omega_y omega_z]
% therefore Xdot = [omega_x omega_y omega_z alpha_x alpha_y alpha_z]
Xdot = zeros(6,1);
Xdot(1,:) = omega_x;
Xdot(2,:) = omega_y;
Xdot(3,:) = omega_z;
Xdot(4,:) = (omega_y*omega_z*(Iyy-Izz))/Ixx;
Xdot(5,:) = (omega_z*omega_x*(Izz-Ixx))/Iyy;
Xdot(6,:) = (omega_x*omega_y*(Ixx-Iyy))/Izz;
end