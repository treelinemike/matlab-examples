% independence, correlation, etc.
% see vector calc notebook for venn diagram

% restart
close all; clear all; clc;

% settings for Latex
% see: http://undocumentedmatlab.com/blog/getting-default-hg-property-values
% get(0,'factory') % return a list of all 'factory...' handles
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultColorbarTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

% other options
doAnimate = 0;
doSaveFrames = 0;

% range of each variable
x = -10:0.1:10;
y = -10:0.1:10;

% means
mean_x = 1;
mean_y = 2;

% variances
var_x = 4;
var_y = 1;

% CORRELATION
% Remember: independence implies uncorrelated
% BUT uncorrelated DOES NOT IMPLY independent
% *EXCEPT* in the case of Gaussians (need a good reference for this...)
% https://en.wikipedia.org/wiki/Normally_distributed_and_uncorrelated_does_not_imply_independent
corr_xy = .8%.65          ;   % [0,1)
animateCorrs = 0:0.01:0.99;

if(~doAnimate)
    animateCorrs = corr_xy;
end

% generate marginal PDFs
pdf_x = normpdf(x,mean_x,sqrt(var_x));
pdf_y = normpdf(y,mean_y,sqrt(var_y));
[XGRID,YGRID] = meshgrid(x,y);

% joint PDF for independent variables (x,y)
pdf_xy_indep = zeros(size(XGRID));
for xidx = 1:length(x)
    for yidx = 1:length(y)
        pdf_xy_indep(yidx,xidx) = pdf_x(xidx)*pdf_y(yidx);
    end
end

% plot results
figure;
set(gcf,'Position',[4.674000e+02 0113 5.648000e+02 0536]);

% marginal PDF of X
axx = subplot(5,5,2:5);
hold on; grid on;
plot(x,pdf_x,'b','LineWidth',1.6);
set(gca,'XTickLabel',{});

% marginal PDF of Y
axy = subplot(5,5,5*(1:4)+1);
hold on; grid on;
plot(-pdf_y,y,'r','LineWidth',1.6);
set(gca,'YTickLabel',{});
set(gca,'YAxisLocation','right');

% plot joint PDF, animate if requested
for corrIdx = 1:length(animateCorrs)
    
    % get next correlation value
    corr_xy = animateCorrs(corrIdx);
    
    % generate covariance matrix
    cov_xy = corr_xy*sqrt(var_x)*sqrt(var_y);
    covar = [var_x cov_xy; cov_xy var_y];
    
    % eigendecomposition of covariance matrix
    % if x and y are correlated, this helps us find a direction in which
    % they can be decoupled (i.e. principal component analysis)
    [S,Lambda] = eig(covar);
    theta = atan2(S(2,1),S(1,1));
    
    % eigenvalues can be computed directly as follows
    % 0.5*( var_x + var_y + sqrt( (var_x + var_y)^2 -4*(1-corr_xy^2)*var_x*var_y ) )
    % 0.5*( var_x + var_y - sqrt( (var_x + var_y)^2 -4*(1-corr_xy^2)*var_x*var_y ) )
    
    % joint PDF of X and Y in general case
    pdf_xy = zeros(size(XGRID));
    for xidx = 1:length(x)
        for yidx = 1:length(y)
            pdf_xy(yidx,xidx) = mvnpdf([x(xidx) y(yidx)],[mean_x mean_y],covar);
        end
    end
    
    % integrate marginal PDFs and joint PDF to confirm that total probability is 1.0
    % if these don't integrate to zero it could be either
    % * domain of x or y is too narrow (gaussian not 'fully contained')
    % * correlation is too high, not enough samples under joint PDF
%     trapz(x,pdf_x)
%     trapz(y,pdf_y)
%     trapz(x,trapz(y,pdf_xy,2))
    
    % joint PDF of X and Y
    axxy = subplot(5,5,[7:10,12:15,17:20,22:25]);
    hold off; axis equal;
    
    % y = beta0 + beta1*x   where beta1 = sqrt(var_y/var_x)
    % to find beta0, force through 2D mean (mean_x,mean_y)
    % ybar = beta0 + beta1*xbar  => beta0 = ybar - beta1*xbar
    % y = ybar - beta1*xbar + beta1*x
    % y = ybar + beta1*(x-xbar)
    
    plot([mean_x mean_x],[min(y) max(y)],'-','Color',[0.75 0.75 0.75],'LineWidth',2);
    hold on; grid on; axis equal;
    plot([min(x) max(x)],[mean_y mean_y],'-','Color',[0.75 0.75 0.75],'LineWidth',2);
    
    plot(x,mean_y+sqrt(var_y/var_x)*(x-mean_x),'-','Color',[0.75 0.75 0.75],'LineWidth',2);
    plot(x,mean_y+sqrt(var_x/var_y)*(mean_x-x),'-','Color',[0.75 0.75 0.75],'LineWidth',2);
    
    plot(x,mean_y+tan(theta)*(x-mean_x),':','Color',[0.5 0.75 0.5],'LineWidth',2);
    plot(x,mean_y+1/(tan(theta))*(mean_x-x),':','Color',[0.5 0.75 0.5],'LineWidth',2);
    
    % plot joint PDF
%     colormap(flipud(brewermap(60,'RdYlBu')));
    contour(XGRID,YGRID,pdf_xy,logspace(-3,-0.3,10),'LineWidth',1);
   
    % add colorbar
    h = colorbar;
    set(h,'Location','south');
    caxis([0 0.5]);
    
    % plot settings
    set(gca,'XAxisLocation','top');
    h1 = title(sprintf('\\textbf{Correlation: %+0.2f}',corr_xy),'FontWeight','Bold','FontSize',12);
    set(h1,'Position',[0, -11.5],'VerticalAlignment','bottom','HorizontalAlignment', 'center');
%     axis equal;
    set(gca,'XLim',[-10 10]);
    set(gca,'YLim',[-10 10]);
    linkaxes([axx axxy],'x');
    linkaxes([axy axxy],'y');
    
    % save figure if requested
    if(doSaveFrames)
        set(gcf,'units','centimeters')
        set(gcf,'papersize',[8,8])
        set(gcf,'paperposition',[0,0,8,8])      
        saveas(gcf,sprintf('frame%03d.png',corrIdx));
    end
    
end

% reset default interpreter
set(0,'defaultTextInterpreter','tex');
set(0,'defaultAxesTickLabelInterpreter','tex');
set(0,'defaultColorbarTickLabelInterpreter','tex');
set(0,'defaultLegendInterpreter','tex');