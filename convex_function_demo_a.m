% restart
close all; clear all; clc;

% options
convex = 1;
doMakeVideo = 0;

% domains for plotting and inspection
x = -10:0.1:10;
x_inspect = -10:2:10;

% evaluate function and set other parameters
if(convex)
    y = x.^2;
    y_inspect = x_inspect.^2;
    plotTitle = '\bfConvex Function';
    outputFilename = 'convex.mp4';
else
    y = 0.024*x.^4-1.6958*x.^2 + 30;
    y_inspect = 0.024*x_inspect.^4-1.6958*x_inspect.^2 + 30;
    plotTitle = '\bfNon-Convex Function';
    outputFilename = 'nonconvex.mp4';
end

% initialize figure
figure;
hold on; grid on;
plot(x,y,'LineWidth',1.6,'Color',[0 0 0.8]);
plot(x_inspect,y_inspect,'.','MarkerSize',20,'Color',[0 0 0.8]);
ph = plot([nan,nan],[nan,nan],'-','LineWidth',3,'Color',[0 0.8 0]);
title(plotTitle);
xlim([-10 10]);
ylim([0 110]);

% check each chord of function
saveFrameIdx = 0;
for startIdx = 1:length(x_inspect)
    for endIdx = (startIdx+1):length(x_inspect)
        startPt = [x_inspect(startIdx) y_inspect(startIdx)];
        endPt = [x_inspect(endIdx) y_inspect(endIdx)];
        ph.XData = [startPt(1),endPt(1)];
        ph.YData = [startPt(2),endPt(2)];
        
        % check to see whether all points in curve lie below line
        x_line = startPt(1):0.01:endPt(1);
        y_line = startPt(2)+(x_line-x_line(1))*((endPt(2)-startPt(2))/(endPt(1)-startPt(1)));       
        y_interp = interp1(x,y,x_line);
        
        % if any points in function are above chord
        % change chord color to red
        if( nnz(y_interp > y_line) > 0 )
            ph.Color = [0.8 0.0 0.0];
        else
            ph.Color = [0.0 0.8 0.0];
        end
        
        % update figure
        drawnow;
        
        % pause or write frames for video
        if(doMakeVideo)
            saveFrameIdx = saveFrameIdx + 1;
            thisImgFile = sprintf('frame%03d.png',saveFrameIdx);
            saveas(gcf,thisImgFile);
            system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
        else
            pause(0.5);
        end
    end
end

% generate movie with ffmpeg
if(doMakeVideo)
    system(['ffmpeg -y -r 4 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' outputFilename]);
    system('del frame*.png');
end