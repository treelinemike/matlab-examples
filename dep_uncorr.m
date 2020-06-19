% example of uncorrelated but NOT independent (i.e. dependent variables that are uncorrelated)

% restart
close all; clear all; clc;

% settings for Latex
% see: http://undocumentedmatlab.com/blog/getting-default-hg-property-values
% get(0,'factory') % return a list of all 'factory...' handles
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultColorbarTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

% parameters
mean_x = 1;
N = 1000;
a = -1;
b = 1;
freqFactor = 5;

% generate samples of X and Y
X = unifrnd(a,b,N,1);
Y = cos(2*pi*freqFactor*X);

% compute sample covariance and correlation
A = [X Y];
A1 = A - repmat(mean(A,1),size(A,1),1);
covXY = (1/(N-1))*A1'*A1
corrXY = covXY(2,1)/sqrt(covXY(1,1)*covXY(2,2))
r2 = corrXY^2;

% variance of X should be (b-a)^2/12 (uniform distribution)
xvar_pop = (b-a)^2/12
xvar_samp = covXY(1,1)

% plot results
figure;
set(gcf,'Position',[0550 0056 1034 0930]);

axx2 = subplot(5,5,2:5);
histogram(X,-1:0.2:1,'FaceColor','b','Normalization','Probability');
set(gca,'XTickMode','manual');
set(gca,'XTick',[]);
title('\boldmath$X \sim $ \textbf{unif}\boldmath$\left( -1,1 \right)$');

axy2 = subplot(5,5,1+(1:4)*5);
histogram(Y,-1:0.2:1,'FaceColor','r','Normalization','Probability');
set(gca,'XTickMode','manual');
set(gca,'XTick',[]);
xlabel(['\boldmath$Y = \cos(2\pi \cdot ' num2str(freqFactor) ' \cdot X)$']);

axxy2 = subplot(5,5,[7:10,12:15,17:20,22:25]);
plot(X,Y,'.','MarkerSize',20,'Color',[0.6 0 0.6])
% xlabel('\textbf{X}');
% ylabel('\textbf{Y}');
h1 = title(sprintf('\\textbf{Correlation: %+0.3f}',corrXY),'FontWeight','Bold','FontSize',12);
set(gca,'XAxisLocation','top');
set(h1,'Position',[0, -1.15],'VerticalAlignment','bottom','HorizontalAlignment', 'center')

% linkaxes([axx2 axxy2],'x');
% linkaxes([axy2 axxy2],'y');
set(gca,'XLim',[-1 1]);
set(gca,'YLim',[-1 1]);
% axis equal;
grid on;
subplot(5,5,1+(1:4)*5)
set(gca,'view',[-90 90]);

% reset default interpreter
set(0,'defaultTextInterpreter','tex');
set(0,'defaultAxesTickLabelInterpreter','tex');
set(0,'defaultColorbarTickLabelInterpreter','tex');
set(0,'defaultLegendInterpreter','tex');