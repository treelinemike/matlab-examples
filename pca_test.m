% Test three related methods of performing basic PCA
% Author: M. Kokko
% Updated: 23-Aug-2021

% restart
close all; clear all; clc; rng('default');

% generate a random data matrix
M = 8;
N = 7000;
X = randn(M,N);

% compute covariance matrices
C = eye(M) - (1/M)*ones(M);
A = C*X;
COV1 = (1/(M-1))*((A')*A); % PARENTHESES VERY IMPORTANT! 
COV2 = (1/(M-1))*(A*(A')); % PARENTHESES VERY IMPORTANT!

% method #1: eigendecomposition of (A')*A
tic
[V1,D1] = eig(COV1);
[d1,ind1] = sort(diag(D1),'descend');
D1s = D1(ind1,ind1);
V1s = V1(:,ind1);
V1s = V1s*diag(sign(V1s(1,:))); % choose first component positive
toc;

% method #2: eigendecomposition of A*(A')
tic
[V2,D2] = eig(COV2);
[d2,ind2] = sort(diag(D2),'descend');
D2s = D2(ind2,ind2);
V2s = unitvec((A')*V2(:,ind2)); % normalize for unit vector columns
V2s = V2s*diag(sign(V2s(1,:))); % choose first component positive
toc;

% method #3: decompose A using the SVD
tic
[U,S,V] = svd(A);
[d3,ind3] = sort(diag(S),'descend');
D3s = (1/(M-1))*S(ind3,ind3).^2;
V3s = V(:,ind3);
V3s = V3s*diag(sign(V3s(1,:))); % choose first component positive
toc;

% compare eigenvalues and eigenvectors
query_idx = 1;
[D1s(query_idx,query_idx) D2s(query_idx,query_idx) D3s(query_idx,query_idx)]
[V1s(1:8,query_idx) V2s(1:8,query_idx) V3s(1:8,query_idx)]
 