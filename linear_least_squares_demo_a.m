% restart
close all; clear; clc;

% options
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
x_opt_mldivide = A\f

% use MATLAB fitlm()
mdl = fitlm(t,f)

% or do it the long way
x_opt_analytic = (A'*A)\(A'*f)     % x_hat = inv(A'*A)*A'*f
resid_analytic = f-A*x_opt_analytic;
SSE_analytic = resid_analytic'*resid_analytic;
rsq_analytic = (SST-SSE_analytic)/SST


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
    for residIdx = 1:length(resid)
        ph(3) = plot(t(residIdx)*ones(1,2),f(residIdx)-[0 resid(residIdx)],'-','LineWidth',1.6,'Color',[0.8 0 0]);
    end
    
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
    ph(2) = plot(xvals,yvals,'-','LineWidth',3,'Color',[0 0.8 0]);
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
    system('del frame*.png');
end

%% check calculations that mean of data lies on regression line
aa = 1/( (length(t))*sum(t.^2)-(sum(t))^2 ) * [length(t) -sum(t);-sum(t) sum(t.^2) ];
(1/length(t))*[sum(t) length(t)]*aa*[sum(t.*f); sum(f)];

%% compute actual cost function
x1_vals = -0.5:0.1:2;
x2_vals = -5:0.2:5;
[X1,X2] = meshgrid(x1_vals,x2_vals);
J = zeros(size(X1));
for elIdx = 1:numel(X1)
   x = [X1(elIdx) X2(elIdx)]';
   J(elIdx) = calcSSE(x,A,f);
end

figure;
set(gcf,'Position',[3.674000e+02 0217 1.023200e+03 4.608000e+02]);
subplot(2,2,[1,3]);
hold on; grid on;
view([234,44]);
sh = surf(X1,X2,J);
ph2(2) = plot3(x_opt_analytic(1),x_opt_analytic(2),SSE_analytic,'.','MarkerSize',40,'Color',[0.8 0 0.8]);
xlabel('\bfSlope');
ylabel('\bfIntercept');
zlabel('\bfCost');
set(sh,'FaceAlpha',0.6);
set(gca,'DataAspectRatio',[1 4 400])

%% now try optimization



% wrap error function to pass parameters/data
optfunc = @(x)calcSSE(x,A,f);

% optimization over both slope and intercept parameters
[x_opt_fminsearch,sse_opt_fminsearch] = fminsearch(optfunc,[0 0]');


% using fmincon with equality constraint forcing regression
% through mean of data
% note: x0 = [0 0]' does NOT satisfy constraints so fmincon adjusts it...
% or we can start with, say, the constraint minimum norm solution
Aeq = [mean(t) 1];
beq = mean(f);
% x0_con = Aeq'*inv(Aeq*Aeq')*beq;
x0_con = [0 0]';
% opts = optimoptions('fmincon','Display','iter-detailed','PlotFcn','optimplotfval','OutputFcn',@myoutfun);



xc2_vals = (beq-x1_vals*Aeq(1))/Aeq(2);
xc = [x1_vals; xc2_vals];
Jc = zeros(1,size(xc,2));
for cIdx = 1:size(xc,2)
   Jc(cIdx) = calcSSE(xc(:,cIdx),A,f); 
end
ph2(5) = plot3(xc(1,:),xc(2,:),Jc,'-','LineWidth',4,'Color',[0.8 0 0.8]);

opts = optimoptions('fmincon','Display','iter','OutputFcn',@myoutfun);
[x_opt_fmincon,sse_opt_fmincon,~,output] = fmincon(optfunc,x0_con,[],[],Aeq,beq,[],[],[],opts);
allData_con = myoutfun([],[],[],1);
allData_con(1,:) = [];
allData_con(end,:) = [];
disp(allData_con);
ph2(4) = plot3(allData_con(:,1),allData_con(:,2),allData_con(:,3),'.-','LineWidth',3,'Color',[0 0.8 0],'MarkerSize',30);
xlim([ min(x1_vals), max(x1_vals)]);
ylim([ min(x2_vals), max(x2_vals)]);

% using fminunc
% opts = optimoptions('fminunc','Display','iter','PlotFcn','optimplotfval');
opts = optimoptions('fminunc','Display','iter','OutputFcn',@myoutfun);
[x_opt_fminunc,sse_opt_fminunc] = fminunc(optfunc,[0 0]',opts);
allData_unc = myoutfun([],[],[],1);
allData_unc(1,:) = [];
allData_unc(end,:) = [];
disp(allData_unc);
ph2(3) = plot3(allData_unc(:,1),allData_unc(:,2),allData_unc(:,3),'.-','LineWidth',3,'Color',[0.8 0 0],'MarkerSize',30);
ph2(1) = plot3(0,0,calcSSE([0;0],A,f),'.','MarkerSize',50,'Color',[0 0.5 0]);
lh = legend(ph2,'Start','Truth','fminunc','fmincon','Constrained Cost');
set(lh,'Position',[0.3465 0.6915 0.0873 0.1866]);


% plot constrained view
subplot(2,2,2);
hold on; grid on;
xlabel('\bfPosition Along Constraint');
ylabel('\bfCost');
title('\bfConstrained View');
posAlongConstraint = sqrt((xc(1,:)-xc(1,1)).^2 +(xc(2,:)-xc(2,1)).^2 );
optPosAlongConstraint = sqrt((allData_con(:,1)-xc(1,1)).^2 +(allData_con(:,2)-xc(2,1)).^2 )
ph3(2) = plot(posAlongConstraint,Jc,'LineWidth',3,'Color',[0.8 0 0.8]);
plot(norm(x_opt_analytic-xc(:,1)),SSE_analytic,'.','MarkerSize',50,'Color',[0.8 0 0.8]);
ph3(1) = plot(optPosAlongConstraint,allData_con(:,3),'.-','LineWidth',2,'MarkerSize',30,'Color',[0 0.8 0]);
xlim([0,max(posAlongConstraint)]);
legend(ph3,'fmincon','Constrained Cost','Location','NorthWest');
% plot cost vs. iterations
subplot(2,2,4);
semilogy(1:size(allData_unc,1),allData_unc(:,3),'.-','MarkerSize',20,'LineWidth',1.6,'Color',[0.8 0 0]);
hold on; grid on;
semilogy(1:size(allData_con,1),allData_con(:,3),'.-','MarkerSize',20,'LineWidth',1.6,'Color',[0 0.8 0]);
legend('fminunc','fmincon');
xlabel('\bfNumber of Iterations');
ylabel('\bfCost');

% function get details of optimization steps
% TODO: this is really a hack, but it mostly works...
function stop = myoutfun(x,optimValues,state,getData)
    persistent optData
    
    if(nargin == 4 && getData)
        stop = optData;
        optData = [];
    else
        optData(end+1,:) = [x' optimValues.fval];
        stop = false;
    end
end

function sse = calcSSE(x,A,b)
    sse = (b-A*x)'*(b-A*x);
end