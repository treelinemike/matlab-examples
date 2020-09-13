% Demonstrate tool tip calibration using data from Medtronic blunt neuro probe and NDI Polaris Optical Tracking system
%
% Takes a series of 3D positions (x,y,z) and orientation quaternions (q0,q1,q2,q3)
% and finds the best (min SSE) choice of an invarient point in the local coordinate
% system (point about which tool is rotated, i.e. the tip).
%
% Optimization is performed numerically (rather than analytically) via
% fminsearch()
%
% Modified: 20200912
% Author:   Mike Kokko

% restart
close all; clear; clc;

% options
filename  = 'tipcal_demo_data.csv';
doMakeVideo = 0;
outputFilename = 'tipcal.mp4';
len_triad = 50;
x0 = [0 0 0]';   % initial position for optimization, in LOCAL coordinates

% read tracker data from file
fileIn = fopen(filename,'r');

% discard one header line
readResult = fgetl(fileIn);

% extract lines from datafile, neglecting rows with missing data
allData = zeros(0,7);
readResult = fgetl(fileIn);
while ( readResult ~= -1)
    if(~contains(readResult,'MISSING'))
        allData(end+1,:) = cell2mat(textscan(readResult,'%*f %*f %*c %f %f %f %f %f %f %f %*f','Delimiter',',','CollectOutput',1));
    end
    readResult = fgetl(fileIn);
end

% use three manually-selected points from the dataset to illustrate
% optimization
sampPts = [170 189 218];
allData = allData(sampPts,:);

% close file
fclose(fileIn);

% function to pass data to cost function
f = @(x_local)probeError(allData(:,5:7)',allData(:,1:4)',x_local);

% run optimization and extract intermediate results at each iteration
opts = optimset('OutputFcn',@myoutfun);
[xOpt, sse] = fminsearch(f,[0 0 0]',opts);
rmse = sqrt(sse/size(allData,1));
optData = myoutfun([],[],[],1);
optData(1,:) = [];
optData(end,:) = [];

% print results to stdout
fprintf('Probe Tip in Local CS: < %+8.3f, %+8.3f, %+8.3f >\n',xOpt);
fprintf('RMSE (mm):  %8.3f\nSSE (mm^2): %8.3f\n',rmse,sse);

% configure plot
figure;
hold on; grid on;
for pointIdx = 1:size(allData,1)
    x = allData(pointIdx,5:7);
    q = allData(pointIdx,1:4)';
    u = quatrotate(q,xOpt);
    
    % compute and show triad
    vx = len_triad*quatrotate(q,[1 0 0]');
    vy = len_triad*quatrotate(q,[0 1 0]');
    vz = len_triad*quatrotate(q,[0 0 1]');    
    plot3([x(1) x(1)+vx(1)],[x(2) x(2)+vx(2)],[x(3) x(3)+vx(3)],'-','LineWidth',3,'Color',[0.8 0 0]);
    plot3([x(1) x(1)+vy(1)],[x(2) x(2)+vy(2)],[x(3) x(3)+vy(3)],'-','LineWidth',3,'Color',[0 0.8 0]);
    plot3([x(1) x(1)+vz(1)],[x(2) x(2)+vz(2)],[x(3) x(3)+vz(3)],'-','LineWidth',3,'Color',[0 0 0.8]);

    % show origin
    plot3(x(1),x(2),x(3),'k.','MarkerSize',30);

    % show optimal point
    plot3(x(1)+u(1),x(2)+u(2),x(3)+u(3),'.','MarkerSize',30,'Color',[0 0.8 0]);
    
    % movable line showing current estimate
    ph(pointIdx) = plot3(x(1)+[0 x0(1)],x(2)+[0 x0(2)],x(3)+[0 x0(3)],'.-','MarkerSize',20,'LineWidth',1.6,'Color',[0.8 0 0.8]);
end
axis image;
view([-18,34]);
xlim([-50 300]);
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');

% now animate to show optimization trajectory
for optIdx = 1:size(optData,1)
    xOptStep = optData(optIdx,1:3)';
    rmseStep = sqrt(optData(optIdx,4)/size(allData,1));
    for pointIdx = 1:size(allData,1)
        x = allData(pointIdx,5:7);
        q = allData(pointIdx,1:4)';
        u = quatrotate(q,xOptStep);
        ph(pointIdx).XData = x(1)+[0 u(1)];
        ph(pointIdx).YData = x(2)+[0 u(2)];
        ph(pointIdx).ZData = x(3)+[0 u(3)];
    end
    title(sprintf('\\bfStep #%03d, RMSE = %6.2f mm',optIdx,rmseStep));
    drawnow;
    
    % pause or write frames for video
    if(doMakeVideo)
        thisImgFile = sprintf('frame%03d.png',optIdx);
        saveas(gcf,thisImgFile);
        system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
    else
        pause(0.1);
    end
end

% generate movie with ffmpeg
if(doMakeVideo)
    system(['ffmpeg -y -r 15 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' outputFilename]);
%     system('del frame*.png');
end



% emperical evaluation of cost function
% there may be a closed-form solution to this
% trying to find the value of x_local that minimizes J (RMSE)
function J = probeError(allX,allQ,x_local)

% first compute global endpoint positions under assumed
% local endpoint position x0
allPoints = zeros(size(allX));
for obsIdx = 1:size(allX,2)
    x = allX(:,obsIdx);  % extract position
    q = allQ(:,obsIdx);  % extract orientation quaternion
    allPoints(:,obsIdx) = x+quatrotate(q,x_local);
end

% compute sum of squared displacements from mean
p = allPoints - sum(allPoints,2)/size(allPoints,2);
J = trace(p*p');
end

% hack to get data out of optimizer
% function get details of optimization steps
% TODO: this is really a hack, but it mostly works...
function stop = myoutfun(x,optimValues,state,getData)
    persistent optData
    
    if(nargin == 4 && getData)
        stop = optData;
        optData = [];
    else
        optData(end+1,:) = [x' optimValues.fval];
        stop = false;
    end
end


% perform quaternion rotation
% essentially using angle (theta) and axis (u, a unit vector) of rotation
% q =     q0     +     q123
% q = cos(theta) + sin(theta)*u
% note: MATLAB includes a similar function in the aerospace toolbox, but
% this is not part of the Dartmouth site license
function v_out = quatrotate(q_in,v_in)

% extract scalar and vector parts of quaternion
q0   = q_in(1);   % real (scalar) part of quaternion
q123 = q_in(2:4); % imaginary (vector) part of quaternion

% rotate v_in using point rotation: v_out = q * v_in * conj(q)
% q_out = quatmult(q_in,quatmult([0; v_in],quatconj(q_in)))
% v_out = q_out(2:4)
% Simplification from Kuipers2002: Quaternions and Rotation Sequences: A Primer with Applications to Orbits, Aerospace and Virtual Reality
v_out = (q0^2-norm(q123)^2)*v_in + 2*dot(q123,v_in)*q123 + 2*q0*cross(q123,v_in);

end