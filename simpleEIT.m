% restart
close all; clear; clc;

% options
doMakeVideo = 1;
outputFilename = 'simpleEIT.mp4';

% define problem
R_true = [100, 120, 112, 87, 132, 97, 101, 99, 114, 31, 94, 103, 87, 111, 102,99]';

A = zeros(9,16);
A( 1,[ 1, 2, 5, 6]) = 1;
A( 2,[ 2, 3, 6, 7]) = 1;
A( 3,[ 3, 4, 7, 8]) = 1;
A( 4,[ 5, 6, 9,10]) = 1;
A( 5,[ 6, 7,10,11]) = 1;
A( 6,[ 7, 8,11,12]) = 1;
A( 7,[ 9,10,13,14]) = 1;
A( 8,[10,11,14,15]) = 1;
A( 9,[11,12,15,16]) = 1;
A = 0.01*A;

V_meas = [
    4.4903
    4.3055
    4.0010
    3.7554
    3.2309
    3.9551
    3.4226
    3.3694
    4.0035];

% prepare figure
figure;
set(gcf,'Position',[0029 2.314000e+02 1444 0420]);
t = tiledlayout(1,4);
t.Padding = 'compact';  % 'normal', 'compact', or 'none'
t.TileSpacing = 'compact';  % 'normal', 'compact', or 'none'

% compute minimum norm solution
x_hat_minnorm = A'*((A*A')\V_meas);

reg_params = logspace(-6,-0.5,100);
for regIdx = 1:length(reg_params)
    
    % compute L1 regularized solution
    x0 = zeros(size(x_hat_minnorm));
    f_L1 = @(x)L1_objective_a(x,A,V_meas,reg_params(regIdx));
    options = optimoptions('fminunc','Display','off');
    [x_hat_L1, U_opt] = fminunc(f_L1,x0,options);
    
    % comput Tikhonov solution
    x_hat_L2 = inv(A'*A + reg_params(regIdx)*eye(16))*(A')*V_meas;
    
    % plot results
    solns = [R_true,x_hat_minnorm,x_hat_L1,x_hat_L2];
    solnLabels = {'Truth','Min. Norm',sprintf('L1 ( \\mu = %06.2e )',reg_params(regIdx)),sprintf('L2 ( \\mu = %06.2e )',reg_params(regIdx))};
    for solnIdx = 1:size(solns,2)
        ax(solnIdx) = nexttile(solnIdx);%subplot(1,3,solnIdx);
        hold on; grid on;
        % surf(XX,YY,reshape(R_true,4,4));
        imagesc(reshape(solns(:,solnIdx),4,4)');
        colorbar('Location','SouthOutside');
        axis image;
        caxis([40 160]);
        ax(solnIdx).XTickLabel = [];
        ax(solnIdx).YTickLabel = [];
        ax(solnIdx).XColor = 'none';
        ax(solnIdx).YColor = 'none';
        title(['\bf' solnLabels{solnIdx}]);
    end
    
    % save frame
    drawnow;
    
    % pause or write frames for video
    if(doMakeVideo)
        thisImgFile = sprintf('frame%03d.png',regIdx);
        saveas(gcf,thisImgFile);
        system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
    end
    
end

% generate movie with ffmpeg
if(doMakeVideo)
    system(['ffmpeg -y -r 15 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' outputFilename]);
    system('del frame*.png');
end

% objective function for L1 regularization
function J = L1_objective_a(x,A,b,regParam)
J = (A*x-b)'*(A*x-b) + regParam*sum(abs(x));
end

