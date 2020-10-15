% restart
close all; clear all; clc;

% color order for plotting
colors = [0 0 .9; .9 0 0; 0 .65 0; .75 0 .75; 0 .8 .8; .8 .6 0; 1 .3 0];

% figure handles
fig_3D = 1;
fig_camview = 2;

% pinhole camera model parameters
fd = 0.05;   % focal length

% define scene in absolute (world) coordinates
scenePts = 2*[
    0  -1   0   0;
    1   0   0   0;
    0   1   0   0;
    0   0   1   0;
    0   0   0   0 ]';

% define translational component of camera trajectory
% using homogeneous coordinates
y = -10:2:10;
x = 4*ones(size(y));
z = 2*ones(size(y));
camTrajT = [x;y;z;zeros(size(y))];
camTrajQ = zeros(4,size(camTrajT,2));

% find rotation & quaterntion taking global xyz to camera xyz
% assume the camera -z axis is directed toward origin at all times
% and orient camera such that camera x axis is parallel to the global
% xy plane
u = -1*camTrajT(1:3,:) ./ vecnorm(camTrajT(1:3,:));
camTrajR = zeros(9,size(u,2));
for camPosIdx = 1:size(camTrajT,2)
    thisU = u(:,camPosIdx);
    thisCamPos = camTrajT(1:3,camPosIdx);
    
    camZ = thisU;
    
    % cross orientation unit vector with global z to get camera x, then
    % rescale to be a unit vector and make sure it has the correct sign
    camX = cross(thisU,[0 0 1]');
    camX = camX./vecnorm(camX);
    
    % finally cross camera z into camera x to get camera y
    camY = cross(camZ,camX);
    camY = camY./vecnorm(camY);
    
    % compile rotation matrix
    % convert to quaternion
    R = [camX camY camZ];
    camTrajR(:,camPosIdx) = reshape(R,9,1);
    
    % TODO: convert rotation matrix to quaternion
    q = rot2quat(R);
    camTrajQ(:,camPosIdx) = q;
end

% display scene
figure(fig_3D);
set(gcf,'Position',[0325 2.994000e+02 0560 4.216000e+02]);
hold on; grid on;
for featureIdx = 1:size(scenePts,2)
    thisColor = colors(mod(featureIdx-1,size(colors,1))+1,:);
    plot3(scenePts(1,featureIdx),scenePts(2,featureIdx),scenePts(3,featureIdx),'.','MarkerSize',20,'Color',thisColor);
end

% display camera path
plot3(x,y,z,'r-','LineWidth',3);
axis equal;
    
% display position, orientation, and camera image plane view for each point
% of interst along trajectory
for camPosIdx = 1:size(camTrajT,2)   
    
    % load camera position and orientation from trajectory
    thisCamPos = camTrajT(:,camPosIdx);
    thisCamR = reshape( camTrajR(:,camPosIdx),3,3);
    thisCamQ = camTrajQ(:,camPosIdx);

    % compute image plane by simulating perspective projection
    featureImageCoords = getFeatureImageCoords(thisCamPos,thisCamQ,fd,scenePts);
       
    % show camera position
    figure(fig_3D);    
    plot3(thisCamPos(1),thisCamPos(2),thisCamPos(3),'r.','MarkerSize',30);
    plot3(thisCamPos(1)+[0,thisU(1)],thisCamPos(2)+[0,thisU(2)],thisCamPos(3)+[0,thisU(3)],'b-');
    
    % show camera orientation
    figure(fig_3D);
    frameScale = 2;
    plot3(thisCamPos(1) + [0,frameScale*thisCamR(1,1)], thisCamPos(2) + [0,frameScale*thisCamR(2,1)], thisCamPos(3) + [0,frameScale*thisCamR(3,1)],'r--','LineWidth',2);
    plot3(thisCamPos(1) + [0,frameScale*thisCamR(1,2)], thisCamPos(2) + [0,frameScale*thisCamR(2,2)], thisCamPos(3) + [0,frameScale*thisCamR(3,2)],'g--','LineWidth',2);
    plot3(thisCamPos(1) + [0,frameScale*thisCamR(1,3)], thisCamPos(2) + [0,frameScale*thisCamR(2,3)], thisCamPos(3) + [0,frameScale*thisCamR(3,3)],'b--','LineWidth',2);
    view([30.3400,35.1200]);
    
    % show camera view
    figure(fig_camview);
    set(gcf,'Position',[8.954000e+02 2.978000e+02 0560 4.200000e+02]);
    hold off; 
    for featureIdx = 1:size(featureImageCoords,2)
        thisColor = colors(mod(featureIdx-1,size(colors,1))+1,:);
        plot(featureImageCoords(1,featureIdx),featureImageCoords(2,featureIdx),'.','MarkerSize',30,'Color',thisColor);
        set(gca,'YDir','reverse');
        hold on;
    end
    grid on;
    axis image;
    xlim([-.03 .03]);
    ylim([-.03 .03]);
    
    % delay before moving to next position
    pause(0.05);
end
