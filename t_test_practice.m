% check to see if we reject with alpha = 0.05 we actually reject H0
% incorrectly 5% of the time....
n = 10;
mu = 80;
sigma = 20;
results = [];
pvals = [];
for trial = 1:10000
   
    % generate sample
    s = randn(n,1)*sigma+mu;
    xbar = mean(s);
    sx = std(s);
    
    % run two-tailed t test
    t = (xbar-mu) / (sx/sqrt(n));
    if(xbar < mu)
        p = 2*tcdf(t,n-1);
    else
        p = 2*tcdf(t,n-1,'upper');
    end
    
    pvals(end+1) = p;
    
    if p < 0.05
        results(end+1) = 1;
    else
        results(end+1) = 0;
    end
end

sum(results)/length(results)