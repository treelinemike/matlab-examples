% reference http://cs229.stanford.edu/notes/cs229-notes1.pdf

% restart
close all; clear figures; clc;

% logistic kernel
l = @(thisArg) 1./(1+exp(-thisArg));

% number of training examples
N = 50;

% generate data
% each training example is a scalar
% labels are scalars in [0,1]
x = linspace(1,10,N);
y = ~[zeros(1,floor(N/2)) ones(1,N-floor(N/2))];

% add some noise to labels
% flips a random number of labels in random positions
shuffle1 = randperm(N);
noiseIdx = shuffle1(1:binornd(floor(N/10),0.5));
y(noiseIdx) = ~y(noiseIdx);

% initial guess for parameters
theta = [0 0]';

% learning rate
alpha = 0.2;

% initialize new figure
figure;

% run stochastic gradient ASCENT (logistic regression uses max likelihood)
shuffle_order = randperm(N);
x_shuf = x(shuffle_order);
y_shuf = y(shuffle_order);
for i = 1:500   % note: should actually terminate based on convergence (likelihood function value, step size, etc.)
    
    for j = 1:N
        thisX = [1 x_shuf(j)]';
        thisY = y_shuf(j);
        theta = theta + alpha*(thisY - l(theta'*thisX) )*thisX;
    end
    
    
    % compute and plot fitted function with current parameter values
    x_full = linspace(min(x),max(x),10000);
    y_hat = l(theta'*[ones(size(x_full));x_full]);
    
    hold off;
    plot(x,y,'b.','MarkerSize',25)
    hold on; grid on;
    plot(x_full, y_hat,'r-');
    drawnow
    pause(0.01);
    
end

% to predict the class of a new unseen example x_query
% compute l( theta'*x_query), then threshold the output at (say) 0.5
% could produce a ROC curve by varying the threshold and computing
% sensitivity and 1-specificity
theta
threshVal = 0.5;
[~,threshIdx] = min( abs( y_hat - threshVal ));
plot(x_full(threshIdx*ones(1,2)),[0 1],'m-','LineWidth',2);