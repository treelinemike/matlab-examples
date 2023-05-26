% GPS CONVERSIONS EXAMPLE

% NOTE: Data from a GPS module may come as
% DDMM.MMMM (lat) or DDDMM.MMMM (lon) which
% requires conversion to decimal degrees:
% lat_deg = floor(tab.Var1/100);
% lat_min = (tab.Var1 - 100*lat_deg);
% lat = lat_deg + lat_min/60;
% lon_deg = floor(tab.Var2/100);
% lon_min = (tab.Var2 - 100*lon_deg);
% lon = -1*(lon_deg + lon_min/60);

% restart
close all; clear; clc;

% general parameters / assumptions
Re = 3963.19;   % [mi] Earth radius (NOTE UNITS!)
alt = 143;      % [m] altitude above ellipsoid (NOTE UNITS!)

% Occom Pond
q1.lat = 43.7111093;    % [deg] 
q1.lon = -72.2877465;   % [deg]

% Lebanon, NH
q2.lat = 43.6423;       % [deg]
q2.lon = -72.2518;      % [deg]

% % Greenwich, UK
% q1.lat = (51+29/60); % [deg]
% q1.lon = 0; % [deg]
% 
% % Boston, MA USA
% q2.lat = (42+21/60); % [deg]
% q2.lon = -1*(71+4/60); % [deg]

% % on prime meridian
% q1.lat = -5; % [deg]
% q1.lon = 0; % [deg]
% 
% % directly across globe
% q2.lat = 5; % [deg]
% q2.lon = 180; % [deg]

%% 1. Compute distance IN CARTESIAN ECEF SPACE with WGS84 ellipsoid approximation
% note: will be shorter than geodesic as line in cartesian space goes
% through the surface of the earth...
xyz_1 = latlon_to_xyz_ecef(q1.lat,q1.lon,alt);
xyz_2 = latlon_to_xyz_ecef(q2.lat,q2.lon,alt);
dist_wgs84_ecef = norm(xyz_2-xyz_1)*0.000621371 % [miles]

%% 2. Compute distance IN CARTESIAN SPACE on a perfect sphere
[x1,y1,z1] = sph2cart(q1.lon*pi/180,q1.lat*pi/180,Re);
[x2,y2,z2] = sph2cart(q2.lon*pi/180,q2.lat*pi/180,Re);
dist_sph2cart = norm([x1,y1,z1]-[x2,y2,z2])  % [miles]

%% 3. Compute distance along geodesic (shortest great circle arc) using SPHERICAL TRIG
% from Wells "New Plane and Spherical Trigonometry" (1902!)
% distance along surface of Earth (spherical trig)
% compute known sides and included angle
b = (pi/2) - q2.lat*(pi/180);
c = (pi/2) - q1.lat*(pi/180);
A = abs(q2.lon - q1.lon)*(pi/180);

% compute two unknown angles: B, C
K1 = 2*atan2(cot(0.5*A)*sin(0.5*(b-c)),sin(0.5*(b+c)));
K2 = 2*atan2(cot(0.5*A)*cos(0.5*(b-c)),cos(0.5*(b+c)));
B = (0.5*(K1+K2));
C = (0.5*(K2-K1));

% compute one unknown side: a
a = 2*atan2(sin(0.5*(B+C))*tan(0.5*(b-c)),sin(0.5*(B-C)));
dist_geodesic = Re*a % [miles]'
xyz_geodesic = [0 dist_geodesic*cos((pi/2)-B);0 dist_geodesic*sin((pi/2)-B);0 0];  % [miles]

%% 4. Show ECEF & Geodesic results "on a map" in ENU (East/North/Up) coordinates
originlat = q1.lat;
originlon = q1.lon;
Rlat = [1 0 0; 0 cosd(90-originlat) -sind(90-originlat); 0 sind(90-originlat) cosd(90-originlat)];
Rlon = [cosd(90+originlon) -sind(90+originlon) 0; sind(90+originlon) cosd(90+originlon) 0; 0 0 1];
R = (Rlon*Rlat)';
xyz_ecef = ([xyz_1 xyz_2]-xyz_1)*0.000621371; % [miles];
xyz_enu = ((R)*(xyz_ecef));
fprintf("Map distance %0.3f mi\n",norm(xyz_ecef(:,2)));
[d,m,s] = deg2dms(B*180/pi);
fprintf("Bearing q1 to q2: %+04d°%02d'%04.1f''\n",d,m,s);
[d,m,s] = deg2dms(-C*180/pi);
fprintf("Bearing q2 to q1: %+04d°%02d'%04.1f''\n",d,m,s);

figure;
hold on; grid on; axis equal;
plot3(xyz_enu(1,:),xyz_enu(2,:),xyz_enu(3,:),'.-','MarkerSize',20,'LineWidth',2,'Color',[0.8 0 0.8]);
plot3(xyz_geodesic(1,:),xyz_geodesic(2,:),xyz_geodesic(3,:),':','LineWidth',4,'Color',[0 0.8 0]);
xlabel('\bfEast');
ylabel('\bfNorth');
zlabel('\bfUp');
legend('ENU','Geodesic');


%% conversion from LAT/LON to ECEF backed out of code from https://github.com/Stellacore/peridetic
function xyz = latlon_to_xyz_ecef(lat_deg,lon_deg,alt_m)
% params = [6378137.0, 298.257222100883]; % GRS80
ellipsoid_params = [6378137.0, 298.257223563];    % WGS84

aa = ellipsoid_params(1);
ff = 1/ellipsoid_params(2);
bb = (1-ff)*aa;
shape = [aa,bb];

up_at_lpa = [cosd(lat_deg)*cosd(lon_deg); cosd(lat_deg)*sind(lon_deg); sind(lat_deg)];

% original shape
radA = shape(1);
radB = shape(2);
lambda = sqrt(radA*radB);
musqs = [(radA)^2; (radA)^2; (radB)^2];

% normalized shape
normPerOrig = 1/lambda;
radA_n = normPerOrig*radA;
radB_n = normPerOrig*radB;
lambda_n = sqrt(radA_n*radB_n);
musqs_n = [(radA_n)^2; (radA_n)^2; (radB_n)^2];

sum_up_mu_sq = dot(up_at_lpa.^2,musqs_n);
scl = lambda/sqrt(sum_up_mu_sq);

xyz = [ (scl*musqs_n(1) + alt_m) * up_at_lpa(1); ...
    (scl*musqs_n(2) + alt_m) * up_at_lpa(2); ...
    (scl*musqs_n(3) + alt_m) * up_at_lpa(3) ];

end

%% conversion from angle to degrees/minutes/seconds
% acutally exists in mapping toolbox... degrees2dms
function [d,m,s] = deg2dms(ang_deg)
    deg_rem = rem(ang_deg,1);
    d = ang_deg-deg_rem;
    ang_min = 60*abs(deg_rem);
    min_rem = rem(ang_min,1);
    m = ang_min-min_rem;
    s = 60*min_rem;
end