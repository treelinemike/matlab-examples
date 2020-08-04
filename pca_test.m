% restart
close all; clear all; clc;

% define test matrix
% pts_ex1 = [1 1 1; 2 1 3; 4.5 1 1; 7 -2 4; 6 3 -2; 4 8 -2; -10 2 5];
% pts_ex2 = [0 -2 7; 3 1 1.3; 8 4 0.3; 4 2.3 11; 7 0.1 -18; 6.5 7 -11; -2.8 4 -9];
% pts_ex3 = [3 3 21; 4 -1 5; -2 3.34 17; 1 16 -2; 0.9 3 4.2; 15 0.2 12; -2.1 -6 5.6];
% pts_ex4 = [5 -5 4.1; 5 -1 9; 6 1 -2.5; 9 -8 -2; 1 11.2 6; -12 3.6 14; 7.4 7 2];
% 
% x = [ pts_ex1(:), pts_ex2(:), pts_ex3(:), pts_ex4(:) ]; 

x = [1 2 3 4; 3 1 6 4; 5 -2 4 -2]';
x = randn(7000,8);

s = size(x,2);
C = eye(s) - (1/s)*ones(s);
A = x*C;
COV1 = (1/(s-1))*(A*(A')); % NOTE: PARENTHESES ARE VERY IMPORTANT AROUND x*x' !!! WILL GO COMPLEX WITHOUT THEM!
COV2 = (1/(s-1))*((A')*A);

% run eigendecomposition

tic
[V1,D1] = eig(COV1);
[d1,ind1] = sort(diag(D1),'descend');
D1s = D1(ind1,ind1);
V1s = V1(:,ind1);
V1s = V1s*diag(sign(V1s(1,:))); % preference for first component of eigenvectors to be positive
toc;

tic
[V2,D2] = eig(COV2);
[d2,ind2] = sort(diag(D2),'descend');
D2s = D2(ind2,ind2);
V2s = unitvec(A*V2(:,ind2));
V2s = V2s*diag(sign(V2s(1,:)));
toc;

tic
[U,S,V] = svd(A);
[d3,ind3] = sort(diag(S),'descend');
D3s = S(ind3,ind3);
V3s = unitvec(U(:,ind2));
V3s = V3s*diag(sign(V3s(1,:)));
toc;