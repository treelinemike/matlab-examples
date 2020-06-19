% restart
close all; clear all; clc;

% define parameter sets
Nvals = [1 4 16 64 128 4096];
pvals = [0.01 0.25 0.50 0.75];

% initialize figure
figure;
set(gcf,'Position',[0001 0049 1920 0946]);

% iterate through each combination of parameters
% generating each PMF
for Nidx = 1:length(Nvals)
    for pidx = 1:length(pvals)
        i = i+1;
        N = Nvals(Nidx);
        p = pvals(pidx);
        pmass = zeros(1,N+1);
        xvals = 0:N;
        for xidx = 1:length(xvals)
            x = xvals(xidx);
            thisP = (factorial(N)/(factorial(x)*factorial(N-x))) *(p^x)*(1-p)^(N-x);
            pmass(xidx) = thisP;
        end
        
        % ensure PMF sums to 1
        if( abs(sum(pmass)-1) > 10^-3)
            warning(sprintf('bin(%d,%0.2f) does not sum to 1!',N,p));
        end
        
        % plot PMF for this combination
        subplot(length(pvals),length(Nvals),length(Nvals)*(pidx-1)+Nidx);
        hold on; grid on;
        stem(xvals,pmass,'b.','MarkerSize',5);
        plot( N*p*ones(1,2),get(gca,'YLim'),'r-');
        xlabel('\bfx');
        ylabel('\bfP');
        title(sprintf('\\bfN = %d, P = %0.2f',N,p));
        set(gca,'XLim',[0 N]);
        set(gca,'YLim',[0,max(pmass)]);
    end
end
