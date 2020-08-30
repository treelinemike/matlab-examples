% restart
close all; clear; clc;

% choose a starting point
x0 = 4.6;

% define problem
charge_locs = [0, 5];
charge_vals = [2, 10]';
q_test = 1;            
k = 1;                        % this chould be Coulomb's constant 8.99eE9

% define state space
xvals = 0:0.01:5;   

% compute ground truth cost
U = (k*q_test*charge_vals(1)./(xvals)) + (k*q_test*charge_vals(2)./(charge_locs(2)-xvals));

% evaulate cost at initial conditions
y0 = charge_PE_1D(x0,charge_locs,charge_vals,q_test,k);

% allow passing additional parameters to fminsearch objective function
f = @(x)charge_PE_1D(x,charge_locs,charge_vals,q_test,k);

% run unconstrained optimization
[x_opt, U_opt] = fminsearch(f,x0);
[x_opt_unc, U_opt_unc] = fminunc(f,x0);

% plot results
figure;
semilogy(xvals,U,'-','LineWidth',1.6,'Color',[0.0 0.0 0.8]);
hold on; grid on;
ph = semilogy(x0,y0,'.','MarkerSize',30,'Color',[0 0.8 0]);
ph(end+1) = plot(x_opt,U_opt,'.','MarkerSize',60,'Color',[0.8 0 0]);
ph(end+1) = plot(x_opt_unc,U_opt_unc,'.','MarkerSize',30,'Color',[0 0.8 0.8]);
legend(ph,{'Start','fminsearch','fminunc'});
title('\bf1D Coulomb Repulsion Energy Minimization');
xlabel('\bfPosition');
ylabel('\bfPotential Energy');

% objective function for optimization
function test_charge_potential = charge_PE_1D(x,charge_locs,charge_vals,q_test,k)
test_charge_potential = (k*q_test*charge_vals(1)/abs(x)) + (k*q_test*charge_vals(2)/abs(charge_locs(2)-x));
end