% COMPARE FORWARD KINEMATICS COMPUTATION FROM DH PARAMETERS
% NOTE: MATLAB Robotics System Toolbox uses standard formulation of 
% DH parameters, not the Khalil 1986 revision.
%
% T_tip_to_base = T_thislink_to_base * (RZ*TZ*TX*RX) * T_prevlink_to_thislink; 
% where "prevlink" is more distal (closer to end effector)

% restart
close all; clear; clc;

% position of base frame
% NOTE: MATLAB doesn't (yet?) support positioning base anywhere other
% than at the origin
% see: https://www.mathworks.com/matlabcentral/answers/696625-change-base-of-rigidbodytree-like-in-seriallink
t0 = [0 0 0]';

% dh paramter table
% NOTE: setFixedTransform ignores joint column (per documentation)
% https://www.mathworks.com/help/robotics/ref/rigidbodyjoint.setfixedtransform.html
% cols: a, alpha (rad), d, theta (rad)
dh =  [0         0         0    0.7854;
    1.4000         0         0    1.5708;
    1.4000         0         0   -0.7854;
    0.2500         0         0         0];

% separate out theta values for joint at sensor zero point (jv_home)
% from joint values taken from sensors indicating displacement from home position (jv_fromhome)
jv_home = dh(:,4);

jv_fromhome =[
   -0.8169;
    0.1504;
    1.1901;
   -1.5708];

jv_total = jv_home + jv_fromhome;

%% Approach 1: MATLAB Robotics System Toolbox framework

% initialize robot
robot = rigidBodyTree;

% add rigid bodies one-by-one
% could easily do this in a loop
body1 = rigidBody('body1');
jnt1 = rigidBodyJoint('jnt1','revolute');
% jnt1.HomePosition = dh(1,4);
setFixedTransform(jnt1,dh(1,:),'dh');
body1.Joint = jnt1;
addBody(robot,body1,'base');

body2 = rigidBody('body2');
jnt2 = rigidBodyJoint('jnt2','revolute');
% jnt2.HomePosition = dh(2,4);
setFixedTransform(jnt2,dh(2,:),'dh');
body2.Joint = jnt2;
addBody(robot,body2,'body1');

body3 = rigidBody('body3');
jnt3 = rigidBodyJoint('jnt3','revolute');
% jnt3.HomePosition = dh(3,4);
setFixedTransform(jnt3,dh(3,:),'dh');
body3.Joint = jnt3;
addBody(robot,body3,'body2');

body4 = rigidBody('body4');
jnt4 = rigidBodyJoint('jnt4','revolute');
% jnt4.HomePosition = dh(4,4);
setFixedTransform(jnt4,dh(4,:),'dh');
body4.Joint = jnt4;
addBody(robot,body4,'body3');

% generate a configuration (i.e. apply joint values)
myConfig = [];
for jnt_idx = 1:robot.NumBodies
    myConfig(jnt_idx).JointName = sprintf('jnt%d',jnt_idx);
    myConfig(jnt_idx).JointPosition = jv_total(jnt_idx);
end


%% Approach 2: Direct transform computation

% updated DH (assuming all revolute)
dh_update = dh;
dh_update(:,4) = dh_update(:,4) + jv_fromhome;

% compute all transforms
TF_links = [];
TF = eye(4);
TF(1:3,4) = t0;
for jnt_idx = 1:size(dh_update,1)
   TF = compute_transform(TF,dh_update(jnt_idx,:));
   TF_links(:,:,jnt_idx) = TF;
end

%% COMPARE RESULTS
TF_matlab = getTransform(robot,myConfig,'body4','base')
TF_direct = TF_links(:,:,4)



function TF_out = compute_transform(TF_in, dh_row)

% initialize output transformation to input transformation
TF_out = TF_in;

% extract values for this particular joint
a           = dh_row(1);
alpha       = dh_row(2);
D           = dh_row(3);
theta       = dh_row(4);

% assemble transformation matrix
RZ = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
RX = [1 0 0 0; 0 cos(alpha) -sin(alpha) 0; 0 sin(alpha) cos(alpha) 0; 0 0 0 1];
TZ = [1 0 0 0; 0 1 0 0; 0 0 1 D; 0 0 0 1];
TX = [1 0 0 a; 0 1 0 0; 0 0 1 0; 0 0 0 1];
T_KHALIL = RX*TX*RZ*TZ;
T_DH = RZ*TZ*TX*RX;   % this is the standard DH convention, but ISI uses modification from Khalil 1986 paper
TF_out = TF_out*T_DH;

end