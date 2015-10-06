function [llow,lmed,lhigh] = diststats95(mydist)

%get middle 95% and median of a distribution
numl = length(mydist);

llow = floor(numl*0.025);
lhigh = floor(numl*0.975);
lsort = sort(mydist,'ascend');
llow = lsort(llow);
lhigh = lsort(lhigh);
lmed = median(mydist);