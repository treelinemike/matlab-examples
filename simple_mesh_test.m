% restart
close all; clear; clc;

% options
doMakeVideo = 0;
outputFilename = 'simple_mesh.mp4';

% define vertices
V = [
    0.0,   0.0,   0.0;
    2.5,   0.0,   0.0;
    0.7,   1.9,   0.0;
    2.9,   1.7,   1.6;
    ];

% define mesh connectivity
F = [
    1 3 2;
    1 2 4;
    2 3 4;
    1 4 3;
    ];

% construct the mesh
thisMesh = triangulation(F,V);

% compute volume enclosed by mesh
meshVol = makmesh_meshvol(thisMesh)

% compute area of each triangle, and total surface area
meshAreas = makmesh_triarea(thisMesh);
meshSA = sum(meshAreas)

A = 0;
faceIdx = 1;
vertexList = [F(faceIdx,:), F(faceIdx,1)];
for vertexIdx = 1:size(F,2)
    vertexList = [F(faceIdx,:), F(faceIdx,1)];
    v1 = V( vertexList(vertexIdx),:);
    v2 = V( vertexList(vertexIdx+1),:);
    A = A + (1/2)*( v1(1)*v2(2) - v1(2)*v2(1) );
end

% display mesh as a patch object
figure;
hold on; grid on;
patch('Vertices',thisMesh.Points,'Faces',thisMesh.ConnectivityList,...
    'EdgeColor','k','FaceColor','flat',...
    'FaceVertexCData',[0.8 0.2 0.2; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5], ...
    'FaceAlpha',0.2);%,'CDataMapping','direct');
plot3(V(:,1),V(:,2),V(:,3),'.','MarkerSize',30,'Color',[0.8 0 0]);
axis equal;
view([-19,34]);
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');


%% check / compare centers

% center of inscribed circle
ictrs = incenter(thisMesh);

% centroid computed via barycentric coordinates
% (scaled in MATLAB? need the 1/3...)
ctrds = barycentricToCartesian(thisMesh,(1:size(F,1))',(1/3)*ones(size(F,1),3));

% show on each face
for faceIdx = 1:size(F,1)

    % centroid could alternatively be computed by averaging vertices...
    ctrd2 = mean(V(F(faceIdx,:),:));

    % show results
    plot3(ictrs(:,1),ictrs(:,2),ictrs(:,3),'.','MarkerSize',30,'Color',[0.8 0 0]);
    plot3(ctrds(:,1),ctrds(:,2),ctrds(:,3),'.','MarkerSize',30,'Color',[0 0.8 0]);
    plot3(ctrd2(:,1),ctrd2(:,2),ctrd2(:,3),'o','MarkerSize',15,'Color',[0 0 0],'LineWidth',2); 
end

%% animate rotating mesh
figure;
hold on; grid on;
axis equal;
colors = [ 0.5 0.5 0.5; 0.8 0 0; 0 0.8 0; 0 0 0.8;];
for faceIdx = 1:size(F,1)
   patch('Vertices',thisMesh.Points,'Faces',thisMesh.ConnectivityList(faceIdx,:),...
    'EdgeColor','k','FaceColor',colors(faceIdx,:), ...
    'FaceAlpha',.7);%,'CDataMapping','direct');
    
   faceVertices = V(F(faceIdx,:),:);
   centroid = mean(faceVertices);
   plot3(centroid(1),centroid(2),centroid(3),'.','MarkerSize',20,'Color',[ 0 0 0]);
   
   % compute face normal
   v1 = faceVertices(2,:)-faceVertices(1,:);
   v2 = faceVertices(3,:)-faceVertices(1,:);
   n = 0.5*unitvec(cross(v1,v2));
   plot3(centroid(1)+[0 n(1)],centroid(2)+[0 n(2)],centroid(3)+[0 n(3)],'-','LineWidth',2,'Color',[0 0.6 0]);
end
for az = 0:359
    view([az 30]);
    drawnow;
    
    % pause or write frames for video
    if(doMakeVideo)
        thisImgFile = sprintf('frame%03d.png',az+1);
        saveas(gcf,thisImgFile);
        system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
    else
        pause(0.01);
    end
end

% generate movie with ffmpeg
if(doMakeVideo)
    system(['ffmpeg -y -r 40 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' outputFilename]);
    system('del frame*.png');
end