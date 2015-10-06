function [r,stderr] = fitdata_pois_data(x,y)

[r,~,stats] = glmfit(x,y,'poisson');
stderr = stats.se(2);