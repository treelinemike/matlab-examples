% restart
close all; clearvars; clc;

% set wheel radius and angular velocity
r = 5;
omega = -0.25;

% compute wheel outline
theta = 0:0.01:2*pi;
x = r*cos(theta);
y = r*sin(theta)+r;

% display wheel outline
figure;
set(gcf,'Position',[4.882000e+02 4.322000e+02 4.184000e+02 3.296000e+02]);
hold on; grid on;
plot(x,y,'-','LineWidth',1.6,'color',[0 0 1]);
plot(0,r,'o','MarkerSize',12,'LineWidth',2,'color',[0 0 1]);

% data storage
queryPts = [];
velVec = [];
accVec = [];

% polar description of query points
thetaRimR = (-pi/2):(pi/6):pi/2;
thetaRimR(end-1) = [];
rimR = sqrt((r.*cos(thetaRimR)).^2+(r.*sin(thetaRimR)+r).^2);

% convert query points to cartesian
% compute velocity vector at each point
% and plot lines of constant velocity
for rp = rimR
   
   phi = asin((rp^2)/(2*r*rp));
   x = rp*cos(phi);
   y = rp*sin(phi);
   queryPts = [queryPts; x,y;-x,y];
   
   % compute velocity vector at this point
   velVec = [velVec; -omega*y omega*x; -omega*y -omega*x];

   % plot line of constant velocity
   if(x ~= 0)
      phiVals = phi:0.01:(pi-phi);
      plot(rp*cos(phiVals),rp*sin(phiVals),'--','color',[0 0.7 0],'LIneWidth',1.6);
   end
   
end

% show vector field
queryPts = [queryPts; [zeros(length(rimR),1) rimR']];
velVec = [velVec; -omega*rimR' zeros(length(rimR),1)];
quiver(queryPts(:,1),queryPts(:,2),velVec(:,1),velVec(:,2),1.125,'color',[0 0.7 0],'LineWidth',1.6);

% show query points
for pointIdx = 1:size(queryPts,1)
   plot(queryPts(pointIdx,1),queryPts(pointIdx,2),'b.','MarkerSize',25); 
end

% adjust plot settings
axis image;
title('\bfLinear Velocity Map for a Rolling Wheel','FontSize',12)
set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});