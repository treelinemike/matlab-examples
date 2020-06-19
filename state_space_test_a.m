% restart
close all; clear; clc;

k  = 5;
m = 2;

s = 1.23;

A = [0 1; -k/m 0];

inv(s*eye(2)-A)

[m*s/(m*s^2+k) m/(m*s^2+k); -k/(m*s^2+k) m*s/(m*s^2+k)] 

% find eigenvectors of A (by hand and w/ MATLAB eig() function)
x1 = [-i*sqrt(m/k)/sqrt(m/k+1); 1/sqrt(m/k+1)];
x2 = [i*sqrt(m/k)/sqrt(1-m/k); 1/sqrt(1-m/k)];

[vec,val] = eig(A);

S= [x1 x2]   % == vec
L = diag([i*sqrt(k/m) -i*sqrt(k/m)]);  % == val

A_star = inv(S)*A*S; % == Lambda (eigenvalue matrix)

trace(A)
trace(A_star)
det(A)
det(A_star)
eig(A)


