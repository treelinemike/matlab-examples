% restart
close all; clear; clc;

D = [0 2; 3 0];
a = [];
b = [];
c = [];
d = [];
for theta = (0:0.01:90)*pi/180

    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    
    D2 = R*D*R'
    a(end+1) = D2(1,1);
    b(end+1) = D2(1,2);
    c(end+1) = D2(2,1);
    d(end+1) = D2(2,2);

end
figure;
hold on; grid on;
plot(a,'r-');
plot(b,'b-');
plot(c,'g--');
plot(d,'m-');
legend('a','b','c','d');


%% another test
syms a b c d;
theta = pi/4;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
T = [0 b; c 0];
T2 = R*T*R'