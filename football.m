% NOTE: CONSIDER REMOVING LATEX INTERPRETER IN LEGENDS B/C MAC AND SOME
% OLDER VERSIONS OF MATLAB THROW ERRORS WITH THAT

% restart
close all; clear all; clc;

% general options
doMakeVideo = 0; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'football';
videoFrameRate = 20; % [frames/sec]

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 0.353;         % [s] simulation end time
dt = 0.005;     % [s] timestep size

% initial conditions (state vector: [phi phi_dot psi psi_dot]')
phi_0     = 0;      % [rad]
phi_dot_0 = 17.78;  % [rad/s]
theta = 45*pi/180;  % [rad] note: theta is constant in torque-free motion
psi_0 = 0;          % [rad]
psi_dot_0 = 22.35;  % [rad/s]
X0 = [phi_0 phi_dot_0 psi_0 psi_dot_0]'; % [rad rad/s rad rad/s]'
X = X0;

% parameters
m = 0.450;  % [kg] mass
kx = 0.05;  % [m]  radius of gyration (x and y)
kz = 0.03;  % [m]  radius of gyration (z)

% assemble inertia tensor
Ixx = m*kx^2;
Iyy = Ixx;
Izz = m*kz^2;
Icm = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];

% data storage
time = [t0];
data = [X0];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % propigate state
    [T,X] = ode45(@(t,X) footballStateProp(t,X,Icm),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

%% Plot time series results
figure;
subplot(2,1,1);
hold on; grid on;
plot(time,data(1,:),'b','LineWidth',1.6);
plot(time,data(3,:),'r','LineWidth',1.6);
legend('$\bf\phi$','$\bf\psi$','FontSize',12,'Interpreter','latex');
xlabel('Time [sec]');
ylabel('Angle [rad]');

subplot(2,1,2);
hold on; grid on;
plot(time,data(2,:),'b','LineWidth',1.6);
plot(time,data(4,:),'r','LineWidth',1.6);
legend('$\bf\dot\phi$','$\bf\dot\psi$','FontSize',12,'Interpreter','latex');
xlabel('Time [sec]');
ylabel('Angular Velocity [rad/s]');
ylim([-30 30]);

%% Develop and display animation of motion

% initialize figure
figure;
hold on; grid on;
firstrun = 1;

% define Cartesian frames
triad_XYZ = 3*[0 0 0; 1 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 1];
triad_xyz_template = .8*triad_XYZ;

% create football
N = 16; % ball "smoothness" parameter; must be a multiple of four!
assert(mod(N,4) == 0,'N must be a multiple of four!');
[ballx,bally,ballz] = sphere(N);             % generate sphere
ball = surf2patch(ballx,bally,ballz,ballz);  % convert to patch object
ball.vertices(:,3) = 1.6*ball.vertices(:,3); % make football oblong
ballvertices_template = ball.vertices;       % save ball template
ballcolor = repmat([repmat([0.5686    0.4039    0.1647],2*N,1); repmat([0.7686    0.6353    0.4471],2*N,1)],N/4,1);  % set ball coloring

% animate each frame of results
saveFrameIdx = 0;
for tIdx = 1:size(data,2)
    
    % extract state at current timestep
    phi = data(1,tIdx);
    phi_dot = data(2,tIdx);
    psi = data(3,tIdx);
    psi_dot = data(4,tIdx);
    
    % compute angular velocity vector (use ihat, jhat, khat components)
    % must be a 3x1 column vector
    omega_xyz = [0; phi_dot*sin(theta); (phi_dot*cos(theta) + psi_dot)];
    Hcm_xyz = Icm*omega_xyz;
    beta = atan(omega_xyz(2)/omega_xyz(3))*180/pi
    
    % compute Euler angle rotations
    R1 = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];          % rotate by psi about z
    R2 = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];  % rotate ball and xyz frame by theta about x
    R3 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];          % rotate ball and xyz frame by phi about Z
    
    % perform rotations
    ballvertices = ballvertices_template*R1*R2*R3;
    triad_xyz = triad_xyz_template*R2*R3;
    omega_XYZ = R3'*R2'*omega_xyz;
    Hcm_XYZ = R3'*R2'*Hcm_xyz;
    
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
    
    % plot ball as patch object
    patch('Faces',ball.faces,'Vertices',ballvertices,'FaceColor','flat','EdgeColor','none','LineWidth',1,'FaceVertexCData',ballcolor);
    
    % finish formatting axes
    axis equal;
    view([145,30]);
    xlabel('\bfX');
    ylabel('\bfY');
    zlabel('\bfZ');
    xlim([-3 3]);
    ylim([-3 3]);
    if(firstrun)
        legend(ph,{'Angular Velocity','Angular Momentum'},'Location','southoutside','AutoUpdate','off');
        firstrun = 0;
    end
    drawnow;
    
        % save frames for video if requested
    if(doMakeVideo)
        thisImgFile = sprintf('frame%03d.png',saveFrameIdx);
        saveFrameIdx = saveFrameIdx + 1;
        saveas(gcf,thisImgFile);
        system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
    end
    
end

% generate movie with ffmpeg
if(doMakeVideo)
    system(['ffmpeg -y -r ' num2str(videoFrameRate) ' -start_number 1 -i frame%03d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' videoFileName '.mp4']);
    system('rm frame*.png');
end

% propagate state
function Xdot = footballStateProp(t,X,Icm)

% recover moments of inertia
I = Icm(1,1);
Izz = Icm(3,3);

% deconstruct state vector
phi     = X(1);
phi_dot = X(2);
psi     = X(3);
psi_dot = X(4);

% construct Xdot from differential equation
% note:     X    = [phi phi_dot psi psi_dot]
% therefore Xdot = [phi_dot phi_ddot psi_dot psi_ddot]
Xdot = zeros(4,1);
Xdot(1,:) = phi_dot;
Xdot(2,:) = 0;
Xdot(3,:) = psi_dot;
Xdot(4,:) = 0;
end