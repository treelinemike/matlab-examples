% compute dual sphere trajectory using direct forcing model
% requires https://www.mathworks.com/matlabcentral/fileexchange/38550-catenary-hanging-rope-between-two-points

% restart
close all; clear; clc;

% restart and make sure Mac systems can find ffmpeg, etc.
close all; clear; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

global Ft alpha;  % this is a hack! but it should work to get tension and alpha out of ODE sovler

% sim options
dt = 0.01;  % [s]
N_catenary = 100;  % # points in catenary
anim_step = 5;                       % skip this many frames to speed up animation
doMakeVideo = 0;                       % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFrameRate = 20;                   % [frames/sec]
videoFileName = 'dual_sphere';

% general setup
params.l_cable = 15; % [m]
params.g = 9.81;     % [m/s^2]

% sphere A
params.ma = 30;   % [kg]
rax0 = 1;  % [m]
ray0 = 8;  % [m]
vax0 = -15;  % [m/s]
vay0 = 15;  % [m/s]

% sphere B
params.mb = 30;   % [kg]
rbx0 = 0;  % [m]
rby0 = 0;  % [m]
vbx0 = 30; % [m/s]
vby0 = 15;  % [m/s]

% initial setttings
X0 = [ rax0 ray0 rbx0 rby0 vax0 vay0 vbx0 vby0 ]';
X = X0;
t = 0;
done_flag = false;
mode = 0;

% data storage
X_data = [X0; 0; 0];  % Ft and alpha are last two rows
mode_data = [mode];

% run simulation
while(~done_flag)

    % calculate timestep for ODE solving
    odeTime = [t t+dt];

    % run ODE solver
    switch(mode)
        case 0
            [T,X] = ode45(@(t,X) state_prop_separate(t,X,params),odeTime,X);
            X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
            Ft = 0;
            alpha = 0;
        case 1
            [T,X] = ode45(@(t,X) state_prop_tethered(t,X,params),odeTime,X);
            X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    end

    % determine whether we need to switch modes or end sim
    if( mode == 0 && (norm([(X(2)-X(4)),(X(1)-X(3))]) >= params.l_cable) )

        % extract from state vector
        rax = X(1);
        ray = X(2);
        rbx = X(3);
        rby = X(4);
        vax = X(5);
        vay = X(6);
        vbx = X(7);
        vby = X(8);

        % look at velocities along and across cable
        vai_cart = [vax vay]';
        vbi_cart = [vbx vby]';
        theta = atan2((ray-rby),(rax-rbx));
        R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        vai_polar = R*vai_cart;
        vbi_polar = R*vbi_cart;

        % conserve momentum for entire system in direction of cable
        % and for each spere individually in direction normal to cable
        % initial state (i) taken instantaneously before cable snaps tight
        % final state (f) taken instantaneously after cable snaps tight
        vf = (params.ma * vai_polar(1) + params.mb * vbi_polar(1)) / (params.ma + params.mb);
        vaf_polar = [vf vai_polar(2)]';
        vbf_polar = [vf vbi_polar(2)]';
        vaf_cart = (R')*vaf_polar;
        vbf_cart = (R')*vbf_polar;

        % update velocities in state vector
        % to reflect effect of cable impulse
        X(5) = vaf_cart(1);
        X(6) = vaf_cart(2);
        X(7) = vbf_cart(1);
        X(8) = vbf_cart(2);
        
        mode = 1;
%         done_flag = true;
    elseif( mode == 1  && (X(2) < 0 || X(4) < 0))
        done_flag = true;

    end

    % save data
    X_data(:,end+1) = [X; Ft; alpha];
    mode_data(end+1) = mode;


end

% initialize figure
figure;
set(gcf,'Position',[0611 0140 0560 0782]);
subplot(5,1,1:3);
hold on; grid on;
axis equal;
xlim([-10 40]);
ylim([-10 40]);
ph_cable = plot(nan(N_catenary,1),nan(N_catenary,1),'LineWidth',1.6,'Color',[0 0 0]);
ph_a = plot(nan,nan,'.','MarkerSize',30,'Color',[0.8 0 0]);
ph_b = plot(nan,nan,'.','MarkerSize',30,'Color',[0 0 0.8]);

subplot(5,1,4);
hold on; grid on;
plot(X_data(9,:),'-','LineWidth',1.6,'Color',[0.8 0.8 0]);
ph_ft = plot(nan,nan,'.','MarkerSize',30,'Color',[0.8 0.8 0]);
ylabel('\bfTension [N]');

subplot(5,1,5);
hold on; grid on;
plot(X_data(10,:),'-','LineWidth',1.6,'Color',[0.8 0 0.8]);
ph_alpha = plot(nan,nan,'.','MarkerSize',30,'Color',[0.8 0 0.8]);
ylabel('\bfAlpha [rad/s^2]');

subplot(5,1,1:3);
% step through data biltting out the new plot positions
saveFrameIdx = 0;
for step_idx = 1:anim_step:size(X_data,2);

    % extract from state vector
    rax = X_data(1,step_idx);
    ray = X_data(2,step_idx);
    rbx = X_data(3,step_idx);
    rby = X_data(4,step_idx);

    % update cable with catenary or straight line depending upon mode
    if(mode_data(step_idx) == 0)
        [x_rope,y_rope] = catenary([rax,ray],[rbx,rby],params.l_cable,N_catenary);
        ph_cable.XData = x_rope;
        ph_cable.YData = y_rope;
    else
        ph_cable.XData = [rax rbx]';
        ph_cable.YData = [ray rby]';
    end

    % update points
    ph_a.XData = rax;
    ph_a.YData = ray;
    ph_b.XData = rbx;
    ph_b.YData = rby;

    plot((params.ma*rax+params.mb*rbx)/(params.ma+params.mb),(params.ma*ray+params.mb*rby)/(params.ma+params.mb),'.','MarkerSize',10,'Color',[0 0.8 0]);

    ph_ft.XData = step_idx;
    ph_ft.YData = X_data(9,step_idx);
    ph_alpha.XData = step_idx;
    ph_alpha.YData = X_data(10,step_idx);

    % update plot
    drawnow;

    % save frame for video if desired
    if(doMakeVideo)
        thisImgFile = sprintf('frame%03d.png',saveFrameIdx);
        saveFrameIdx = saveFrameIdx + 1;
        saveas(gcf,thisImgFile);
        system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
    else
        pause(0.05);
    end
end

% generate movie with ffmpeg if desired
if(doMakeVideo)
    system(['ffmpeg -y -r ' num2str(videoFrameRate) ' -start_number 1 -i frame%03d.png -vf "format=rgba,scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' videoFileName '.mp4']);
    system('rm frame*.png');
end


% state propagation function for tethered sphere motion
function Xdot = state_prop_tethered(t,X,params)

global Ft alpha;

% extract parameters
g = params.g;

% deconstruct necessary components of the current state vector
rax = X(1);
ray = X(2);
rbx = X(3);
rby = X(4);
vax = X(5);
vay = X(6);
vbx = X(7);
vby = X(8);

% generate and solve linear system for soln = [aax, aay, abx, aby, Ft, alpha]'
theta = atan2((ray-rby),(rax-rbx));
omega = (vbx-vax)/(rby-ray);

A = [ ...
    params.ma, 0, 0, 0, cos(theta), 0; ...
    0, params.ma, 0, 0, sin(theta), 0 ; ...
    0, 0, params.mb, 0, -cos(theta), 0; ...
    0, 0, 0, params.mb, -sin(theta), 0; ...
    -1, 0, 1, 0, 0, (ray-rby); ...
    0, -1, 0, 1, 0, (rbx-rax)];
b = [ 0, -params.ma*params.g, 0, -params.mb*params.g, omega*(vby-vay), omega*(vax-vbx)]';
soln = A\b;
Ft = soln(5);
alpha = soln(6);

% construct Xdot from differential equations
Xdot = zeros(8,1);
Xdot(1) = vax;
Xdot(2) = vay;
Xdot(3) = vbx;
Xdot(4) = vby;
Xdot(5) = soln(1);  % aax
Xdot(6) = soln(2);  % aay
Xdot(7) = soln(3);  % abx
Xdot(8) = soln(4);  % aby

end



% state propagation function for separate sphere motion
function Xdot = state_prop_separate(t,X,params)

% extract parameters
g = params.g;

% deconstruct necessary components of the current state vector
vax = X(5);
vay = X(6);
vbx = X(7);
vby = X(8);

% construct Xdot from differential equations
Xdot = zeros(8,1);
Xdot(1) = vax;
Xdot(2) = vay;
Xdot(3) = vbx;
Xdot(4) = vby;
Xdot(5) = 0;   % aax
Xdot(6) = -g;  % aay
Xdot(7) = 0;   % abx
Xdot(8) = -g;  % aby

end