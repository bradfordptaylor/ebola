function [r,stderr] = fitdata_pois(x,y)

%fit data with exponential with poisson distributed noise
[p,~,stats] = glmfit(x,y,'poisson');
r = p(2);   %best fit
stderr = stats.se(2);   %standard error