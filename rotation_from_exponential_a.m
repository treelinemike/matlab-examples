% restart
close all; clear; clc;

% define rotation
theta = 25*pi/180;
u = [0.3842, 0.6533, 0.6524]';

% construct rotation matrix from quaternion
q = [cos(theta/2); sin(theta/2)*u];
Rq = quat2matrix(q)

% construct rotation matrix from exponential (via Taylor series for expm(); see MATLAB doc for expm() which includes references to better? methods)
A = theta*[0 -u(3) u(2);u(3) 0 -u(1);-u(2) u(1) 0];
Re = zeros(3);
errStorage = [];
for n = 0:5
   Re = Re + (1/factorial(n))*A^n;
   err = abs(Rq-Re);
   errStorage(n+1) = max(err(:));
end
Re

% MATLAB version
Rm = expm(A)

plot(errStorage);