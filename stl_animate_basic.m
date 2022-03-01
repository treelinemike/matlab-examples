% Test of animating an STL file to follow a simple trajectory with rigid
% transformations
%
% System requirements (on Windows) for making video w/ ffmpeg
% - ffmpeg: https://www.gyan.dev/ffmpeg/builds/ (get "essentials" release build, extract somewhere, add the ffmpeg/bin path to PATH environmental variable
% - imagemagick: https://imagemagick.org/script/download.php... be sure to install legacy utils for *convert*
% - gnuwin coreutils: http://gnuwin32.sourceforge.net/packages/coreutils.htm.... and add to PATH, alternatively just change system('rm.... to system('del...
%
% Notes: 
% - STL units should be mm
% - For Solidworks export, use option "Do not translate STL output data to positive space"
%
% Author:   Mike Kokko
% Modified: 01-Mar-2022

% restart
close all; clear; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
anim_step = 5;                      % skip this many steps between frames to speed up animation
doMakeVideo = 0;                    % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'stl_animtate';
videoFrameRate = 30;                % [frames/sec]

% load STL file (see notes above)
stl_raw = stlread('candycane.stl');

% reduce mesh size for better display
[f,v] = reducepatch(stl_raw.ConnectivityList,stl_raw.Points,0.1);

% physical parameters taken from SolidWorks
m   = 1.459;        % [kg] mass
Ixx = 0.001002;     % [kg*m^2] moment of inertia about body x axis
Iyy = 0.004685;     % [kg*m^2] moment of inertia about body y axis
Izz = 0.005523;     % [kg*m^2] moment of inertia about body z axis
Icm  = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];  % inertia matrix taken at CM about principal axes
R_pa = [0.44 0.90 0; -0.90 0.44 0; 0 0 1];
x_cm = [0.024340, 0.108379, 0.0];

% rotate STL file to its principal axes and move CM to origin
stl_pts = ((v/1000)-x_cm)*R_pa';
stl = triangulation(f,stl_pts);

% initialize figure
figure;
hold on; grid on; axis equal;
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
view([30 30]);

% define an arbitrary path in space
traj.t = 0:0.01:10;
traj.x = traj.t*(0.25/traj.t(end));
traj.y = zeros(size(traj.t));
traj.z = 0.0625*sin((2*pi/(traj.t(end))) .* traj.t);
traj.grad = gradient(traj.z,traj.x);
plot3(traj.x,traj.y,traj.z,'-','LineWidth',2,'Color',[0 0 0.8]);

% initialize STL plotting
ph = patch('Faces',stl.ConnectivityList,'Vertices',nan(size(stl.Points)),'FaceColor',[0.8 0.2 0.2],'EdgeColor',[0 0 0],'LineWidth',0.5);
saveFrameIdx = 0;

for t_idx = 1:size(traj.t,2)

    % only display if we're a the correct step spacing
    if( mod(t_idx-2,anim_step) == 0 )

        % compute and apply appropriate rotation
        theta = atan(traj.grad(t_idx)) + pi/2;  % [rad]
        R = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
        ph.Vertices = stl.Points*R' + [traj.x(t_idx) traj.y(t_idx) traj.z(t_idx)]; % note: cm shift included

        % update plot
        xlim([-0.1 0.35]);
        ylim([-0.1 0.1]);
        zlim([-0.20 0.15]);        
        drawnow;

         % save frames for video if requested
        if(doMakeVideo)
            thisImgFile = sprintf('frame%04d.png',saveFrameIdx);
            saveFrameIdx = saveFrameIdx + 1;
            saveas(gcf,thisImgFile);
            system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
        end
    end

end

% generate movie with ffmpeg
if(doMakeVideo)
    system(['ffmpeg -y -r ' num2str(videoFrameRate) ' -start_number 1 -i frame%04d.png -vf "format=rgba,scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' videoFileName '.mp4']);
    system('rm frame*.png');
end