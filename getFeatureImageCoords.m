% simulate a camera snapshot of scene feature points with known global
% positions given by (targetLocs)
function featureImageCoords = getFeatureImageCoords(camPos, camQuat, camFL, featureLocs)

featureImageCoords = zeros(size(featureLocs));
R = quat2rot(camQuat);

% If the recorded image was simply a projection of
% the scene points onto the camera xy plane we could use the code below.
% We actually need a PERSPECTIVE projection
% A = [camX camY];
% camViewPts = ((A'*A)\A')*scenePts(1:3,:)

%%%%%% METHOD 1: VECTOR ANALYSIS APPROACH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% camZ = R(:,3);
% for featureIdx = 1:size(featureLocs,2)
%     % get unit vector pointing to this target from the camera origin
%     A = featureLocs(1:3,featureIdx) - camPos(1:3);
%     Auv = A./vecnorm(A);
%     if abs(abs(dot(Auv,camZ))-1) > 0.0001
%         x = (camFL/dot(Auv,camZ))*Auv - camFL*camZ;
%     else
%         x = [0 0 0]';
%     end
%     camViewPoint = R\x;
%     assert( abs(camViewPoint(3)) < 0.00001, 'Camera geometry error.');
%     featureImageCoords(:,featureIdx) = [camViewPoint;0];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% METHOD 2: MATRIX APPROACH STARING FROM POINT IN IMAGE FRAME %%%%%
% REF: Scaramuzza and Fraundorfer 2011 (IEEE Robotics & Automation Mag.)
% TODO: allow u0, v0 to be passed as parameters so image origin can be
%       placed in upper left hand corner as is standard in image processing
u0 = 0;
v0 = 0;
A = featureLocs(1:3,:) - camPos(1:3);
uv = [camFL 0 0; 0 camFL 0; u0 v0 1]*(R\A);
featureImageCoords = uv ./ uv(3,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end