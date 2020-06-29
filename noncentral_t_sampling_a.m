% restart
close all; clear all; clc;

% options
doMakeMovies = 1;
frameCount = 1;  % initialize at 1 for ffmpeg

%% simulate the sampling distribution of the t-statistic
mu = 2;
sigma = 2;
N = 3;

sampDist = [];
statDist = [];

x = -10:0.001:10;

figure;
set(gcf,'Position',[0488 2.426000e+02 4.402000e+02 0372]);

numSamples = 1e5;
for sampIdx = 1:numSamples;
    
    thisSamp = mu + sigma*randn(N,1);
    xbar = mean(thisSamp);
    s = std(thisSamp);
    
    sampDist(end+1,1) = xbar;
    statDist(end+1,1) = (xbar)/(s/sqrt(length(thisSamp)));
    
    if(sampIdx <= 50 || sampIdx == numSamples)
        subplot(3,1,1);
        hold off;
        plot(x,normpdf(x,mu,sigma),'b-','LineWidth',1.6);
        hold on; grid on;
        plot(thisSamp,zeros(size(thisSamp)),'.','MarkerSize',25,'Color',[0 0.7 0]);
        xlim([-5,5]);
        % ylim([0 1]);
        title(sprintf('\\bf\\mu = %0.0f    \\sigma = %0.0f    N = %d,  Samples = %6d',mu,sigma,N,sampIdx));
        plot(xbar,0,'r.','MarkerSize',25);
        plot(xbar*ones(2,1),get(gca,'Ylim'),'r-','LineWidth',1.2);
        
        subplot(3,1,2);
        hold off;
        plot(sampDist,zeros(size(sampDist)),'b.','MarkerSize',25);
        hold on; grid on;
        plot(x,ksdensity(sampDist,x),'b-','LineWidth',1.6);
        plot(xbar,0,'r.','MarkerSize',25);
        plot(xbar*ones(2,1),get(gca,'Ylim'),'r-','LineWidth',1.2);
        xlim([-5,5]);
        % ylim([0 1]);
        
        subplot(3,1,3);
        hold off;
        plot(x,normpdf(x,0,1),'--','LineWidth',1.6,'Color',[0 0 .7]);
        hold on; grid on;
        plot(statDist,zeros(size(sampDist)),'m.','MarkerSize',25);
        plot(x,ksdensity(statDist,x),'m-','LineWidth',1.6);
        plot(statDist(end),0,'r.','MarkerSize',25);
        % plot(statDist(end)*ones(2,1),get(gca,'Ylim'),'r-','LineWidth',1.2);
        xlim([-5,5]);
        % ylim([0 1]);
        
        
        
        % save frames for animation, or pause
        % then convert - not in MATLAB although could use system() command - using a codec compatible with LaTeX (Beamer)
        % see: https://tex.stackexchange.com/questions/141622/trouble-using-ffmpeg-to-encode-a-h264-video-in-a-mp4-container-to-use-with-medi
        % first download FFMPEG for windows: https://ffmpeg.zeranoe.com/builds/ OR https://www.ffmpeg.org/download.html
        % ffmpeg -r 10 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 output.mp4
        if(doMakeMovies)
            saveas(gcf,sprintf('frame%03d.png',frameCount));
            frameCount = frameCount + 1;
        else
            pause(0.05);
        end
        
        
    end
end

% save animation
if(doMakeMovies)
    system(['ffmpeg -y -r 4 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 sampling_dist_of_t_statistic.mp4']);
    system('del frame*.png');
end

%% compare t and non-central t
figure

mu = 2;
s = 2;
N = 3;
SE = s/sqrt(N);
theta = mu/SE;

x = -5:0.01:5;
y1 = tpdf(x,9);
y2 = nctpdf(x,9,theta);

hold on;
plot(x,y1,'b-','LineWidth',1.6);
plot(x,y2,'r-','LineWidth',1.6);

%% NONCENTRAL T POWER ANALYSIS

doMakeMovies = 1;
frameCount = 1;  % initialize at 1 for ffmpeg

figure
set(gcf,'Position',[4.082000e+02 3.898000e+02 0560 2.168000e+02]);
delta = 10;     % [mm]
sigma = 11.06;  % [mm]
N = 1;          % initialize to 1 so first value is actually 2 (N-1 must be > 0)
x = -10:0.01:10;
truePower = 0;

while ( truePower < 0.8)
    
    N = N + 1;
    
    theta = abs(delta/(sigma/sqrt(N)));  % non-centrality parameter
    threshold = tinv(0.95,N-1);
    truePower = 1-nctcdf(threshold,N-1,theta);
    
    % plot null hypothesis distribution
    h1 = tpdf(x,N-1);
    h2 = nctpdf(x,N-1,theta);
    
    
    hold off;
    plot(x,h1,'b-','LineWidth',1.6);
    hold on; grid on;
    plot(x,h2,'r-','LineWidth',1.6);
    
    thresholdIdx = find(x >= threshold,1,'first');
    powerX = nctinv(0.2,N-1,theta);
    powerIdx = find(x >= powerX,1,'first');
    area([x(thresholdIdx) x(thresholdIdx:end)],[0 h1(thresholdIdx:end)],'LineStyle','none','FaceColor','b','FaceAlpha',0.2);
    area([x(powerIdx) x(powerIdx:end)],[0 h2(powerIdx:end)],'LineStyle','none','FaceColor','r','FaceAlpha',0.2);
    plot(x(thresholdIdx)*ones(1,2),[0,0.4],'--','LineWidth',1.6,'Color',[0 0.7 0]);
    
    ylim([0 0.4]);
    title(sprintf('\\bfN = %d   Power = %0.4f',N,truePower),'FontSize',12);
    drawnow
    
    % save frames for animation, or pause
    % then convert - not in MATLAB although could use system() command - using a codec compatible with LaTeX (Beamer)
    % see: https://tex.stackexchange.com/questions/141622/trouble-using-ffmpeg-to-encode-a-h264-video-in-a-mp4-container-to-use-with-medi
    % first download FFMPEG for windows: https://ffmpeg.zeranoe.com/builds/ OR https://www.ffmpeg.org/download.html
    % ffmpeg -r 10 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 output.mp4
    % the -r 10 means 10 frames per second (??), try -r 1.2 for this...
    if(doMakeMovies)
        saveas(gcf,sprintf('frame%03d.png',frameCount));
        frameCount = frameCount + 1;
    else
        pause(0.5);
    end
    
end

% save animation
if(doMakeMovies)
    system(['ffmpeg -y -r 4 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 non_central_t_power.mp4']);
    system('del frame*.png');
end

%% **CENTRAL** T POWER ANALYSIS

doMakeMovies = 1;
frameCount = 1;  % initialize at 1 for ffmpeg

figure
set(gcf,'Position',[4.082000e+02 3.898000e+02 0560 2.168000e+02]);
delta = 10;     % [mm]
sigma = 11.06;  % [mm]
N = 1;          % initialize to 1 so first value is actually 2 (N-1 must be > 0)
x = -10:0.01:10;
truePower = 0;

while ( truePower < 0.8)
    
    N = N + 1;
    
    t_delta = delta/(sigma/sqrt(N));
    threshold = tinv(0.95,N-1);
    truePower = 1-tcdf(threshold-t_delta,N-1);
    
    % plot null hypothesis distribution
    h1 = tpdf(x,N-1);
    h2 = tpdf(x-t_delta,N-1);
    
    
    hold off;
    plot(x,h1,'b-','LineWidth',1.6);
    hold on; grid on;
    plot(x,h2,'r-','LineWidth',1.6);
    
    thresholdIdx = find(x >= threshold,1,'first');
    powerX = t_delta + tinv(0.2,N-1);
    powerIdx = find(x >= powerX,1,'first');
    area([x(thresholdIdx) x(thresholdIdx:end)],[0 h1(thresholdIdx:end)],'LineStyle','none','FaceColor','b','FaceAlpha',0.2);
    area([x(powerIdx) x(powerIdx:end)],[0 h2(powerIdx:end)],'LineStyle','none','FaceColor','r','FaceAlpha',0.2);
    plot(x(thresholdIdx)*ones(1,2),[0,0.4],'--','LineWidth',1.6,'Color',[0 0.7 0]);
    
    ylim([0 0.4]);
    title(sprintf('\\bfN = %d   Power = %0.4f',N,truePower),'FontSize',12);
    drawnow
    
    % save frames for animation, or pause
    % then convert - not in MATLAB although could use system() command - using a codec compatible with LaTeX (Beamer)
    % see: https://tex.stackexchange.com/questions/141622/trouble-using-ffmpeg-to-encode-a-h264-video-in-a-mp4-container-to-use-with-medi
    % first download FFMPEG for windows: https://ffmpeg.zeranoe.com/builds/ OR https://www.ffmpeg.org/download.html
    % ffmpeg -r 10 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 output.mp4
    % the -r 10 means 10 frames per second (??), try -r 1.2 for this...
    if(doMakeMovies)
        saveas(gcf,sprintf('frame%03d.png',frameCount));
        frameCount = frameCount + 1;
    else
        pause(0.5);
    end
    
end

% save animation
if(doMakeMovies)
    system(['ffmpeg -y -r 4 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 central_t_power.mp4']);
    system('del frame*.png');
end