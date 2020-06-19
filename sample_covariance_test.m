%% compare various methods for computing the covariance matrix from a data sample

% restart
close all; clear all; clc;

% sample data 
x = [161   203   235   176   201   188   228   211   191   178];
y = [159   206   241   163   197   193   209   189   169   201];

%% element-wise computation (very slow, but very clear)
% first compute sum of squares terms
sx2 = 0;
sy2 = 0;
sxy = 0;
n = length(x);
for idx = 1:n
    sx2 = sx2 + (x(idx)-mean(x))^2;
    sy2 = sy2 + (y(idx)-mean(y))^2;
    sxy = sxy + (x(idx)-mean(x))*(y(idx)-mean(y));
end

%% normalize, using n-1 to make unbiased estimator
sx2 = (1/(n-1))*sx2;
sy2 = (1/(n-1))*sy2;
sxy = (1/(n-1))*sxy;
cov1 = [sx2 sxy; sxy sy2]

%% mean correction method
A = [x' y'];
Ac = A - repmat(mean(A,1),n,1);
cov2 = (1/(n-1))*Ac'*Ac

%% mean correctio method from Ali S. Hadi book
% A = data matrix
% C = centering matrix, essentially removes mean
A = [x' y'];
C = diag((1-1/n)*ones(1,n)) + (-1/n)*ones(n,n).*(~diag(ones(1,n)));
cov3 = (1/(n-1))*A'*C*A

%% built-in MATLAB function
cov4 = cov(x,y)  % or equivalently, A = [x' y']; cov4 = cov(A)