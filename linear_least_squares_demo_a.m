% restart
close all; clear; clc;

% options
doShowWeightedSoln = 1;
doRunAnimation = 0;
doMakeVideo = 0;
outputFilename = 'lsq_demo_a.mp4';

% input data
t = [2 5 6 7]';
f = [3 5 4 6]';

% compute coefficient matrix
A = [t ones(size(t))];

% the regression line will always pass through the mean of the data!
meanpt = mean([t f]);
SST = (f-meanpt(2))'*(f-meanpt(2));  % SST = SSE for model that is just mean of data

% solve with MATLAB mldivide()
x_hat_1 = A\f

% use MATLAB fitlm()
mdl = fitlm(t,f)

% or do it the long way
x_hat_2 = (A'*A)\(A'*f)     % x_hat = inv(A'*A)*A'*f
resid_2 = f-A*x_hat_2;
SSE_2 = resid_2'*resid_2;
rsq_2 = (SST-SSE_2)/SST

% compute the weighted solution
COV_f = (1/(length(f)-1))*(f-mean(f))*(f-mean(f))';
WR = inv(COV_f);  % ideally compute this using eigendecomposition, or similar


% initialize plot
figure;
set(gcf,'Position',[0488 2.962000e+02 5.866000e+02 4.658000e+02]);
hold on; grid on;
axis equal;
xlim([0 10]);
ylim([0 10]);

%% now run a simulation for demonstration
% keep track of errors
min_sse = SST;
min_rsq = 0;   % R^2 at mean model = (SST-SST)/SST = 0
min_xhat = [0 meanpt(2)];

% step through angles, pinning solution at mean of data
if(doRunAnimation)
    theta_vals = [0:0.01:pi/2-0.01,pi/2+0.01:0.01:pi];
else 
%     theta_vals = atan(x_hat_1(1));
    theta_vals = 0.2;
end
for stepIdx = 1:length(theta_vals);
    theta = theta_vals(stepIdx);
    % compute trial slope and intercept
    m = tan(theta);
    b = meanpt(2)-m*meanpt(1);

    % compute x_hat vector
    % and resuduals
    x_hat = [m;b];
    f_hat = A*x_hat;
    resid = f-f_hat;

    % compute SSE and update minimum if this is the best so far
    SSE = resid'*resid;
    if(SSE < min_sse)
        min_sse = SSE;
        min_xhat = x_hat;
        min_rsq = (SST-SSE)/SST;  % guaranteed to be in [0,1] only when x_hat minimizes (A*t-f)'*(A*t-f)
    end

    % plot residuals
    cla;
%     for residIdx = 1:length(resid)
%         ph(3) = plot(t(residIdx)*ones(1,2),f(residIdx)-[0 resid(residIdx)],'-','LineWidth',1.6,'Color',[0.8 0 0]);
%     end
    
    % determine how to plot current and best fit lines
    if(m < 1)
       xvals = [0,10];
       yvals = m*xvals+b;
    else
       yvals = [0,10];
       xvals = (yvals-b)/m;
    end
    if(min_xhat(1) < 1)
       min_xvals = [0,10];
       min_yvals = min_xhat(1)*min_xvals + min_xhat(2);
    else
       min_yvals = [0,10];
       min_xvals = (min_yvals-min_xhat(2))/min_xhat(1);
    end
    
    % plot center point
%   plot(meanpt(1),meanpt(2),'.','MarkerSize',30,'Color',[0 0.8 0]);
    
    % plot data
    ph(1) = plot(t,f,'.','MarkerSize',30,'Color',[0 0 0.8]);

    % plot current and best fit lines
%     ph(2) = plot(xvals,yvals,'-','LineWidth',3,'Color',[0 0.8 0]);
    if(doRunAnimation)
        ph(4) = plot(min_xvals,min_yvals,':','LineWidth',3,'Color',[0.8 0 0.8]);
    end
    
    % finish plot
    if(doRunAnimation)
        title(sprintf('SSE = %06.2e, Best SSE = %06.2e, Best R^2 = %04.2f',SSE,min_sse,min_rsq));
        legend(ph,{'Data','Current Model','Residuals','Best Model'},'Location','Northwest');
    else
        legend(ph,{'Data','Model','Residuals'},'Location','Northwest');
    end
    drawnow;
    
    % pause or write frames for video
    if(doMakeVideo)
        thisImgFile = sprintf('frame%03d.png',stepIdx);
        saveas(gcf,thisImgFile);
        system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
    else
        pause(0.01);
    end
    
end

% generate movie with ffmpeg
if(doMakeVideo)
    system(['ffmpeg -y -r 30 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' outputFilename]);
%     system('del frame*.png');
end

%% check calculations that mean of data lies on regression line
aa = 1/( (length(t))*sum(t.^2)-(sum(t))^2 ) * [length(t) -sum(t);-sum(t) sum(t.^2) ];
(1/length(t))*[sum(t) length(t)]*aa*[sum(t.*f); sum(f)];