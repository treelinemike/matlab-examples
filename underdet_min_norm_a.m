% find the minimum norm solution to Ax = b
% equivalently: minimize norm of x subject to the constraint Ax = b
% use lagrange multipliers

% restart
close all; clear; clc;

% constraint
A = [5 1];
b = 4.5;

x_opt_analytical = (A')*((A*A')\b);


% create meshed grid
x1_vals = -3:0.5:3;
x2_vals = -4:0.5:4;
[X1,X2] = meshgrid(x1_vals,x2_vals);

% evaulate f and g on grid
F = (X1.^2+X2.^2);
G = A(1)*X1+A(2)*X2;

x_parallel_grad = x1_vals([1,end]);
y = x_parallel_grad/5;
y_constraint = x2_vals([1,end]);
x_constraint = (4.5-y_constraint)/5;
line_param = 0:0.01:1;
x_const_cost = (1-line_param)*x_constraint(1) + line_param*x_constraint(end);
y_const_cost = (1-line_param)*y_constraint(1) + line_param*y_constraint(end);
z_const_cost = x_const_cost.^2+y_const_cost.^2;


figure;
set(gcf,'Position',[2.778000e+02 0342 7.702000e+02 0420]);
subplot(1,2,1);
hold on; grid on;
surf(X1,X2,F);
surf(X1,X2,G,'EdgeColor','none','FaceColor',[0.8 0 0],'FaceAlpha',0.5);
plot3(x_constraint,y_constraint,b*ones(1,2),'-','LineWidth',3,'Color',[0.8 0 0]);
plot3(x_const_cost,y_const_cost,z_const_cost,'-','LineWidth',5,'Color',[0.8 0 0.8]);
plot3(x_opt_analytical(1),x_opt_analytical(2),norm(x_opt_analytical)^2,'.','MarkerSize',50,'Color',[0.8 0 0.8]);
xlabel('\bfx_1');
ylabel('\bfx_2');
set(gca,'DataAspectRatio',[1 1 5]);
view([-78 34]);
zlim([-10 30]);
lh(1) = legend('f(x)','g(x)','Constraint','Constrained Cost','Optimal Point');
lh(1).Position = [0.3059    0.7012    0.1887    0.2048];


subplot(1,2,2);
hold on; grid on; axis equal;
contour(X1,X2,F,10,'LineColor',[0 0 0.8]);
contour(X1,X2,G,10,'LineColor',[0.8 0 0]);
plot(x_parallel_grad,y,'-','LineWidth',3,'Color',[0 0.8 0]);
plot(x_constraint,y_constraint,'-','LineWidth',3,'Color',[0.8 0 0]);
plot(x_opt_analytical(1),x_opt_analytical(2),'.','MarkerSize',50,'Color',[0.8 0 0.8]);
xlabel('\bfx_1');
ylabel('\bfx_2');
lh(2) = legend('f(x) Level Curves','g(x) Level Curves','Parallel Gradients','Constraint','Optimal Point');
lh(2).Position = [0.7726    0.7012    0.1908    0.2048];
% lh(2).Position = [0.7726    0.7486    0.1908    0.2048];