% restart
close all; clear all; clc;

xi = unifrnd(0,1,10000,1);
yi = norminv(xi);

figure;
set(gcf,'Position',[2.802000e+02 3.354000e+02 0768 4.264000e+02]);

sp1 = subplot(2,2,1);
hold on; grid on;
[a,b] = hist(yi);
a = (length(b)*(a/sum(a)))/(max(b)-min(b));
barh(b,a);
[pdf,pdf_dom] = ksdensity(yi);
plot(pdf,pdf_dom,'r-','LineWidth',1.6);

sp2 = subplot(2,2,2);
hold on; grid on;
x_inv_cdf = 0:0.001:1;
y_inv_cdf = norminv(x_inv_cdf);
plot(x_inv_cdf,y_inv_cdf,'-','Color',[0 0.7 0],'LineWidth',1.6);

sp4 = subplot(2,2,4);
hold on; grid on;
[a2,b2] = hist(xi);
a2 = length(b2)*a2/sum(a2);
bar(b2,a2);
[pdf2,pdf_dom2] = ksdensity(xi);
plot(pdf_dom2,pdf2,'b-','LineWidth',1.6);

linkaxes([sp1 sp2],'y');
linkaxes([sp2 sp4],'x');
drawnow;