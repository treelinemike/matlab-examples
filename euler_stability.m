% NOTE: CONSIDER REMOVING LATEX INTERPRETER IN LEGENDS B/C MAC AND SOME
% OLDER VERSIONS OF MATLAB THROW ERRORS WITH THAT

% simulate Euler's equations of motion for a rectangular board

% restart
close all; clear all; clc;

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 5;         % [s] simulation end time
dt = 0.05;      % [s] timestep size

% physical parameters of board
rho = 0.4e3;                         % [kg/m^3] density
l_x = 0.05;                          % [m] thickness along x
l_y = 0.25;                          % [m] thickness along y
l_z = 0.01;                          % [m] thickness along z
V   = l_x*l_y*l_z;                   % [m^3] volume
m   = rho*V;                         % [kg] mass
Ixx = (1/12)*m*(l_y^2+l_z^2);        % moment of inertia about body x axis
Iyy = (1/12)*m*(l_x^2+l_z^2);        % moment of inertia about body y axis
Izz = (1/12)*m*(l_x^2+l_y^2);        % moment of inertia about body z axis
Icm  = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];  % inertia matrix taken at CM about principal axes

% initial conditions X0 = [theta_x_0 theta_y_0 theta_z_0 omega_x_0 omega_y_0 omega_z_0]
X0 = [0 0 0 2 .2 0]'; % [rad rad rad rad/s rad/s rad/s]'
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
h = 0.7*max([l_x, l_y, l_z]); % scaling parameter
triad_XYZ = h*[0 0 0; 1 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 1];
triad_xyz_template = .8*triad_XYZ;

% assemble model of board in body coordintes
board.v = [ ...
        l_x/2, -l_y/2,  l_z/2;
        l_x/2,  l_y/2,  l_z/2;
        l_x/2,  l_y/2, -l_z/2;
        l_x/2, -l_y/2, -l_z/2;
       -l_x/2, -l_y/2,  l_z/2;
       -l_x/2,  l_y/2,  l_z/2;
       -l_x/2,  l_y/2, -l_z/2;
       -l_x/2, -l_y/2, -l_z/2 ];
board.f = [ ...
    1, 2, 3, 4, 1;
    4, 3, 7, 8, 4;
    8, 7, 6, 5, 8;
    5, 1, 2, 6, 5;
    5, 1, 4, 8, 5;
    6, 2, 3, 7, 6 ];
board.c = repmat([0.2 0.4 0.6 0.8 0 0]',1,3); % face colors

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
        
    % transform point cloud to correct location in inertial space
    board.vrot = board.v*R';
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
    
    % plot board as patch object
    patch('Faces',board.f,'Vertices',board.vrot,'FaceColor','flat','EdgeColor','none','LineWidth',1,'FaceVertexCData',board.c);
    
    % finish formatting axes
    axis equal;
    xlim([-h h]);
    ylim([-h h]);
    zlim([-h h]);
    xlabel('\bfx');
    ylabel('\bfy');
    zlabel('\bfz');
    
    if(firstrun)
        legend(ph,{'Angular Velocity','Angular Momentum'},'Location','southoutside','AutoUpdate','off');
        firstrun = 0;
    end
    view([145,30]);
    drawnow;
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