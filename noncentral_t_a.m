% restart
close all; clear all; clc;

true_mu = 2;
true_sigma = 1;
N = 3;
num_trials = 1e6;
x = -20:0.001:20;
% for N = 2:2:30
% for true_mu = 0:1:5
samp_dist = zeros(num_trials,1);
for trialIdx = 1:num_trials;
    thisSample = true_mu + true_sigma*randn(N,1);
%     samp_dist_samp_mean(trialIdx) = (mean(thisSample)-true_mu)/(std(thisSample)/sqrt(length(thisSample)));
    samp_dist(trialIdx) = (mean(thisSample))/(std(thisSample)/sqrt(length(thisSample)));

end
subplot(1,2,1);
hold on;
ksdensity(samp_dist,x);
plot(x,nctpdf(x,N-1,true_mu/(true_sigma/sqrt(N))),'r-');    

ylim([ 0 0.4]);
xlim([-8,10]);
subplot(1,2,2);
qqplot(samp_dist);

drawnow;



