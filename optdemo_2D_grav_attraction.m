% restart
close all; clear; clc;

% choose a starting point
x0 = [1.3 0.7]';

% define problem
charge_locs = [ 0.5  0.5;
                0.5  1.0;
                1.0  0.5;
                1.0  1.0 ];
charge_vals = [4, 4, 6, 4]';
q_test = 1;            
k = 1;                        % this chould be Coulomb's constant 8.99eE9

% define state space and mesh it
xvals = 0:0.01:1.5;
yvals = 0:0.01:1.5;
[X,Y] = meshgrid(xvals,yvals);

% compute actual potential function
U = zeros(size(X));
for charge_idx = 1:size(charge_locs,1)
   rvals = sqrt((X-charge_locs(charge_idx,1)).^2 + (Y-charge_locs(charge_idx,2)).^2);
   U = U + k*q_test*charge_vals(charge_idx)./rvals;
end

% allow passing additional parameters to fminsearch objective function
f = @(x)charge_PE_2D(x,charge_locs,charge_vals,q_test,k);

% run unconstrained optimization
[x_opt, U_opt] = fminsearch(f,x0);
[x_opt_unc, U_opt_unc] = fminunc(f,x0);

% solve again with consrained optimization

[x_opt_con, U_opt_con] = fmincon(f,x0,[],[],[0 1],x0(2));

% plot results
figure;
hold on; grid on;
contour(X,Y,U,200);
axis equal;
plot([0 1.5],x0(2)*ones(1,2),'m-','LineWidth',1.6);
ph = plot(x0(1),x0(2),'.','MarkerSize',30,'Color',[0 0.8 0]);
ph(end+1) = plot(x_opt(1),x_opt(2),'.','MarkerSize',60,'Color',[0.8 0 0]);
ph(end+1) = plot(x_opt_unc(1),x_opt_unc(2),'.','MarkerSize',30,'Color',[0 0.8 0.8]);
ph(end+1) = plot(x_opt_con(1),x_opt_con(2),'.','MarkerSize',30,'Color',[0.8 0 0.8]);
legend(ph,{'Start','fminsearch','fminunc','fmincon'});
title('\bf2D Gravitation Energy Minimization');

% objective function for optimization
function test_charge_potential = charge_PE_2D(x,charge_locs,charge_vals,q_test,k)
    test_charge_potential = 0;
    for charge_idx = 1:size(charge_locs,1)
       rval = sqrt( (x(1)-charge_locs(charge_idx,1))^2 + (x(2)-charge_locs(charge_idx,2))^2  );
       test_charge_potential = test_charge_potential + -1*k*q_test*charge_vals(charge_idx)/rval;
    end
end