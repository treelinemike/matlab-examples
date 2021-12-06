% Simulate double sphere penetrometer using Euler's equations of motion & simple ballistic motion 
% Based on ENGS 72 simulations
%
% System requirements (on Windows) for making video w/ ffmpeg
% - ffmpeg: https://www.gyan.dev/ffmpeg/builds/ (get "essentials" release build, extract somewhere, add the ffmpeg/bin path to PATH environmental variable
% - imagemagick: https://imagemagick.org/script/download.php... be sure to install legacy utils for *convert*
% - gnuwin coreutils: http://gnuwin32.sourceforge.net/packages/coreutils.htm.... and add to PATH, alternatively just change system('rm.... to system('del...
%
% Author: Mike Kokko
% Updated: 06-Dec-2021

% restart and make sure Mac systems can find ffmpeg, etc.
close all; clear; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
anim_step = 100;                       % skip this many frames to speed up animation
doMakeVideo = 0;                       % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFrameRate = 20;                   % [frames/sec]
videoFileName = 'dual_sphere';

% simulation time parameters
t0 = 0;                                % [s] simulation start time
tf = 20;                               % [s] simulation end time
dt = 0.001;                            % [s] timestep size

% physical parameters of system
% ASSUMING IDENTICAL SPHERES (TODO: allow different masses)
rho = 7850;                            % [kg/m^3] density
R   = 0.125;                           % [m] sphere radius
r   = 5.00;                            % [m] sphere displacement from CM (half length of cable)
V   = (4/3)*pi*R^3;                    % [m^3] volume
m   = rho*V;                           % [kg] mass of single sphere
Ixx = (4/5)*m*(R^2);                   % moment of inertia about body x axis
Iyy = (4/5)*m*(R^2) + 2*m*(r^2);       % moment of inertia about body y axis
Izz = (4/5)*m*(R^2) + 2*m*(r^2);       % moment of inertia about body z axis
Icm  = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];    % inertia matrix taken at CM about principal axes

% assemble into parameter structure
params.Icm = Icm;
params.g = 1.62;                       % [m/s^2]
params.m1 = m;                         % [kg] assuming equal mass spheres for now...
params.m2 = m;                         % [kg] assuming equal mass spheres for now...
params.sph_locs = [r 0 0; -r 0 0];     % [m] rigid body sphere locations - relative to point "p" on cable

% DEFINE INITAL CONDITIONS
R0 = eye(3);                           % initial rotation matrix
omega0 = [0 -1.8 0];                     % [rad/s] initial rigid body angular velocity
tctr0 = [0 0 0];                       % [m] initial position of point "p" on cable in inertial frame
vctr0 = 15*[cos(pi/4) 0 sin(pi/4)];    % [m/s] initial linear velocity at point "p" on cable

% make sure we're working with the CM, adjust as necessary
% (i.e. if point "p" on cable is not the CM)
cm_offset = (params.m1*params.sph_locs(1,:)+params.m2*params.sph_locs(2,:))/(params.m1+params.m2);
params.sph_locs = params.sph_locs - cm_offset;  % now params.sph_locs is relative to CM!

% compute initial position and velocity at the CM
cm_offset_rot = cm_offset*R0'; 
tcm0 = tctr0 + cm_offset_rot;
vcm0 = vctr0 + cross(omega0,cm_offset_rot);

% run forward kinematics to find ICs at each sphere
% this determines position and velocity of each individual sphere
% consistend with the given initial conditions
t0_rel = params.sph_locs*R0';             % [m] positions of each sphere relative to CM
t1_0 = tcm0 + t0_rel(1,:);                % [m] inital position of first sphere
v1_0 = vcm0 + cross(omega0,t0_rel(1,:));  % [m/s] initial velocity of first sphere
t2_0 = tcm0 + t0_rel(2,:);                % [m] inital position of second sphere
v2_0 = vcm0 + cross(omega0,t0_rel(2,:));  % [m/s] initial velocity of second sphere

% initialize state vector 
X0 = [0, 0, 0, ...                     % theta_x, theta_y, theta_z (rigid body rotation angles, all in local frame! - use diffentially ONLY) 
      omega0(1:3), ...                 % omega_x, omega_y, omega_z (rigid body angular rates - all in local frame!)
      tcm0, ...                        % tcm_x, tcm_y, tcm_z (CM position - in inertial frame!)
      vcm0, ...                        % vcm_x, vcm_y, vcm_z (CM velocity - in inertial frame!)
      t1_0, ...                        % t1_x, t1_y, t1_z (sphere 1 position - in inertial frame!)
      v1_0, ...                        % v1_x, v1_y, v1_z (sphere 1 velocity - in inertial frame!)
      t2_0, ...                        % t2_x, t2_y, t2_z (sphere 2 position - in inertial frame!)
      v2_0, ...                        % v2_x, v2_y, v2_z (sphere 2 velocity - in inertial frame!)
      ]';
X = X0;

% mode indicator
% 1: ballistic flight, connected
% 2: ballistic flight, independent
% could add more...
sim_mode = 1;   

% data storage
time = [t0];                           % [s] track timestamps
data = [X0];                           % track state vector
quat_data = [matrix2quat(R0)];         % track overall rotation matrix compressed in quaternion form
mode_data = [sim_mode];
v1_mag = nan;
v2_mag = nan;

%% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % check for state transition, changing modes if necessary
    cord_uv = unitvec(data(19:21,end) - data(13:15,end));
    P = [1 0 0; 0 1 0; 0 0 0];     % P = A'*A where columns of A are the basis of the plane to project into
    proj_uv = unitvec(P*cord_uv);
    incl = acos(abs(dot(cord_uv,proj_uv)));  % angle of incination above xy plane

    % get height and vertical speed of CM
    tcm_z = data(9,end);
    vcm_z = data(12,end);

    % evaluate switching conditions
    switch(sim_mode)
        case 1
            % switch to individal sphere motion when we are nearly flat and within a set distance of the surface
            if( (incl < (1*pi/180)) && (tcm_z < 30)  && (vcm_z < 0) )
                sim_mode = 2;
            end
        case 2
            % stop each sphere when it contacts regolith
            t1_z = data(15,end);
            t2_z = data(21,end);
            if(t1_z < 0)
                X(16:18) = 0;
                v1_mag = norm( data(16:18) );
                params.m1 = 0;  % hack... we're really want a new mode for motion through regolith
            end
            if(t2_z < 0)
                X(22:24) = 0;
                v2_mag = norm( data(22:24) );
                params.m2 = 0;  % hack... we're really want a new mode for motion through regolith
            end
    end

    % save mode
    mode_data(end+1) = sim_mode;

    % propagate state
    switch(sim_mode)
        case 1   % rigid body motion

            % propagate state
            [T,X] = ode45(@(t,X) state_prop_rb(t,X,params),odeTime,X);
            X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
            
            % compute rotation between body coordinates and inertial space
            theta_local = X(1:3);
            d_theta_x = theta_local(1) - data(1,end);
            d_theta_y = theta_local(2) - data(2,end);
            d_theta_z = theta_local(3) - data(3,end);
            xRot = [1 0 0 ; 0 cos(d_theta_x) -sin(d_theta_x); 0 sin(d_theta_x) cos(d_theta_x)];
            yRot = [cos(d_theta_y) 0 sin(d_theta_y); 0 1 0; -sin(d_theta_y) 0 cos(d_theta_y)];
            zRot = [cos(d_theta_z) -sin(d_theta_z) 0 ; sin(d_theta_z) cos(d_theta_z) 0; 0 0 1];
            R  = quat2matrix(quat_data(:,end))*(xRot*yRot*zRot*eye(3));     % incremental rotations so order shouldn't matter
            quat_data(:,end+1) = matrix2quat(R);

            % use forward kinematics to update the positions and velocities of each of the two spheres
            omega_local = X(4:6)';
            tcm = X(7:9)';
            vcm = X(10:12)';
            t_rel = params.sph_locs*R';                % [m] positions of each sphere relative to CM
            t1 = tcm + t_rel(1,:);                     % [m] position of first sphere
            v1 = vcm + cross(omega_local,t_rel(1,:));  % [m/s] velocity of first sphere
            t2 = tcm + t_rel(2,:);                     % [m] position of second sphere
            v2 = vcm + cross(omega_local,t_rel(2,:));  % [m/s] velocity of second sphere
            X(13:15) = t1;
            X(16:18) = v1;
            X(19:21) = t2;
            X(22:24) = v2;

        case 2 % independent sphere motion
            
            % propagate state
            [T,X] = ode45(@(t,X) state_prop_separate(t,X,params),odeTime,X);
            X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
            X(1:12) = nan;
            quat_data(:,end+1) = quat_data(:,end);
            
    end

    % store results from this timestep
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
    time(end+1) = T(end);

end
fprintf('v1 at impact: %0.2f m/s\n', v1_mag);
fprintf('v2 at impact: %0.2f m/s\n', v2_mag);

%% develop and display animation of motion
% define Cartesian frames
triad_XYZ = 0.8*r*[0 0 0; 1 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 1];
triad_xyz_template = .8*triad_XYZ;

% initialize figure
figure;
set(gcf,'Position',[2.506000e+02 1.914000e+02 1084 4.832000e+02]);
hold on; grid on;

% show full trajectory
plot3(data(7,:),data(8,:),data(9,:),'--','LineWidth',1.0,'Color',[0 0 0.8]);
plot3(data(13,:),data(14,:),data(15,:),'--','LineWidth',1.0,'Color',[0.8 0 0]);
plot3(data(19,:),data(20,:),data(21,:),'--','LineWidth',1.0,'Color',[0.8 0 0]);

% plot cord and spheres
ph.cord = plot3([nan nan],[nan nan],[nan nan],'.-','LineWidth',2,'Color',[0 0 0]);
ph.sph1 = plot3(nan,nan,nan,'o','MarkerSize',10,'MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor','none');
ph.sph2 = plot3(nan,nan,nan,'o','MarkerSize',10,'MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor','none');

% plot triad and grab handles for blitting
ph.XYZ = plot3(triad_XYZ(:,1),triad_XYZ(:,2),triad_XYZ(:,3),'LineWidth',4,'Color','k');
ph.xyz_x = plot3(nan(1,2),nan(1,2),nan(1,2),'LineWidth',3,'Color',[0.7 0 0]);
ph.xyz_y = plot3(nan(1,2),nan(1,2),nan(1,2),'LineWidth',3,'Color',[0 0.7 0]);
ph.xyz_z = plot3(nan(1,2),nan(1,2),nan(1,2),'LineWidth',3,'Color',[0 0 0.7]);

% add additional plot features
th = title('');
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
% view([145,30]);
view([180,0]);

% animate each frame of the results
saveFrameIdx = 0;
for tIdx = 1:size(data,2)
    
    % skip frames if desired to speed up animation
    % don't do this in the for loop b/c need to update rotation at each step
    if( mod(tIdx-1,anim_step) == 0 )

        % extract sphere positions
        % don't blit them out yet for debugging
        t1 = data(13:15,tIdx)';
        t2 = data(19:21,tIdx)';
        t_sph = [t1;t2];

        % mode-specific blitting
        sim_mode = mode_data(tIdx);
        switch(sim_mode)
            case 1
                % get position and rotation
                tcm = data(7:9,tIdx)';
                R   = quat2matrix(quat_data(:,tIdx));

                % blit spheres and cord
                ph.sph1.XData = t_sph(1,1);
                ph.sph1.YData = t_sph(1,2);
                ph.sph1.ZData = t_sph(1,3);
                ph.sph2.XData = t_sph(2,1);
                ph.sph2.YData = t_sph(2,2);
                ph.sph2.ZData = t_sph(2,3);
                ph.cord.XData = t_sph(:,1);
                ph.cord.YData = t_sph(:,2);
                ph.cord.ZData = t_sph(:,3);

                % blit the body xyz frame
                triad_xyz = tcm + triad_xyz_template*R';
                ph.xyz_x.XData = [triad_xyz(1,1) triad_xyz(2,1)];
                ph.xyz_x.YData = [triad_xyz(1,2) triad_xyz(2,2)];
                ph.xyz_x.ZData = [triad_xyz(1,3) triad_xyz(2,3)];
                ph.xyz_y.XData = [triad_xyz(3,1) triad_xyz(4,1)];
                ph.xyz_y.YData = [triad_xyz(3,2) triad_xyz(4,2)];
                ph.xyz_y.ZData = [triad_xyz(3,3) triad_xyz(4,3)];
                ph.xyz_z.XData = [triad_xyz(5,1) triad_xyz(6,1)];
                ph.xyz_z.YData = [triad_xyz(5,2) triad_xyz(6,2)];
                ph.xyz_z.ZData = [triad_xyz(5,3) triad_xyz(6,3)];

            case 2

                % blit spheres and clear cord
                ph.sph1.XData = t_sph(1,1);
                ph.sph1.YData = t_sph(1,2);
                ph.sph1.ZData = t_sph(1,3);
                ph.sph2.XData = t_sph(2,1);
                ph.sph2.YData = t_sph(2,2);
                ph.sph2.ZData = t_sph(2,3);
                ph.cord.XData = nan;
                ph.cord.YData = nan;
                ph.cord.ZData = nan;

                % clear the triad
                ph.xyz_x.XData = nan(1,2);
                ph.xyz_x.YData = nan(1,2);
                ph.xyz_x.ZData = nan(1,2);
                ph.xyz_y.XData = nan(1,2);
                ph.xyz_y.YData = nan(1,2);
                ph.xyz_y.ZData = nan(1,2);
                ph.xyz_z.XData = nan(1,2);
                ph.xyz_z.YData = nan(1,2);
                ph.xyz_z.ZData = nan(1,2);
        end
                               
        % finish formatting axes
        axis equal;
        xlim([-2*r 1.1*max([ data(13,:) data(19,:)])]);
        ylim([-2*r 2*r]);
        zlim([-2*r 1.5*max([ data(15,:) data(21,:)])]);       
        th.String = sprintf('Dual Sphere Sim (%6.3fs)',time(tIdx));

        % save frame for video if desired
        if(doMakeVideo)
            thisImgFile = sprintf('frame%03d.png',saveFrameIdx);
            saveFrameIdx = saveFrameIdx + 1;
            saveas(gcf,thisImgFile);
            system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
        end
        drawnow;
    end
end

% generate movie with ffmpeg if desired
if(doMakeVideo)
    system(['ffmpeg -y -r ' num2str(videoFrameRate) ' -start_number 1 -i frame%03d.png -vf "format=rgba,scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' videoFileName '.mp4']);
    system('rm frame*.png');
end

%% plot trajectories
% could add a bunch of other traces here...

figure;

ax = subplot(3,1,1);
hold on; grid on;
plot(time,data(7,:),'-','Color',[0.8 0 0],'LineWidth',1.6);
plot(time,data(8,:),'-','Color',[0 0.8 0],'LineWidth',1.6);
plot(time,data(9,:),'-','Color',[0 0 0.8],'LineWidth',1.6);
legend('t_x','t_y','t_z');
xlabel('\bfTime [sec]');
ylabel('\bf t [m]');

ax(end+1) = subplot(3,1,2);
hold on; grid on;
plot(time,data(10,:),'-','Color',[0.8 0 0],'LineWidth',1.6);
plot(time,data(11,:),'-','Color',[0 0.8 0],'LineWidth',1.6);
plot(time,data(12,:),'-','Color',[0 0 0.8],'LineWidth',1.6);
legend('v_x','v_y','v_z');
xlabel('\bfTime [sec]');
ylabel('\bf v [m/s]');

ax(end+1) = subplot(3,1,3);
hold on; grid on;
plot(time,data(4,:),'-','Color',[0.8 0 0],'LineWidth',1.6);
plot(time,data(5,:),'-','Color',[0 0.8 0],'LineWidth',1.6);
plot(time,data(6,:),'-','Color',[0 0 0.8],'LineWidth',1.6);
legend('w_x','w_y','w_z');
xlabel('\bfTime [sec]');
ylabel('\bf w [rad/s]');

% % set axis limits
linkaxes(ax,'x');
xlim(time([1 end]));
% ylim([ min([0, min(min(data(4:6,:)))]) max([0, max(max(data(4:6,:)))]) ]);

%% state propagation function for rigid body motion
% only propagates rigid body rotation and translation of CM
% forward kinematics can be used OUTSIDE OF THIS FUNCTION to compute
% positions and velocities of each sphere
function Xdot = state_prop_rb(t,X,params)

% extract parameters 
Icm = params.Icm;
g = params.g;

% recover moments of inertia
Ixx = Icm(1,1);
Iyy = Icm(2,2);
Izz = Icm(3,3);

% deconstruct necessary components of the current state vector
omega_x = X(4);
omega_y = X(5);
omega_z = X(6);
vcm_x   = X(10);
vcm_y   = X(11);
vcm_z   = X(12);

% construct Xdot from differential equation
Xdot = zeros(24,1);
Xdot(1,:) = omega_x;
Xdot(2,:) = omega_y;
Xdot(3,:) = omega_z;
Xdot(4,:) = (omega_y*omega_z*(Iyy-Izz))/Ixx;
Xdot(5,:) = (omega_z*omega_x*(Izz-Ixx))/Iyy;
Xdot(6,:) = (omega_x*omega_y*(Ixx-Iyy))/Izz;
Xdot(7,:) = vcm_x;
Xdot(8,:) = vcm_y;
Xdot(9,:) = vcm_z;
Xdot(10,:) = 0;
Xdot(11,:) = 0;
Xdot(12,:) = -g;
end

%% state propagation function for separate sphere motion
% only propagates rigid body rotation and translation of CM
% forward kinematics can be used OUTSIDE OF THIS FUNCTION to compute
% positions and velocities of each sphere
function Xdot = state_prop_separate(t,X,params)

% extract parameters
g = params.g;

% deconstruct necessary components of the current state vector
v1_x    = X(16);
v1_y    = X(17);
v1_z    = X(18);
v2_x    = X(22);
v2_y    = X(23);
v2_z    = X(24);

% construct Xdot from differential equations
Xdot = zeros(24,1);
Xdot(13,:) = v1_x;
Xdot(14,:) = v1_y;
Xdot(15,:) = v1_z;
Xdot(16,:) = 0;
Xdot(17,:) = 0;
Xdot(18,:) = -g;
Xdot(19,:) = v2_x;
Xdot(20,:) = v2_y;
Xdot(21,:) = v2_z;
Xdot(22,:) = 0;
Xdot(23,:) = 0;
Xdot(24,:) = -g;
end