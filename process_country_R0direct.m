function process_country_R0direct(country,R0max)
%turn trajectories in rate constant distributions and save statistics
 
%country specific number of sampled lambdas as used in clusterfun.m
R0vec = [1.0125 1.025 1.05 1.1:.1:2.5];
R0vecorig = 1.1:.1:2.5;
vecdiff = length(R0vec)-length(R0vecorig);
%loop over lambda values contained in clusterfun.m
lambdadist = cell(1,R0max+1);
zerotimedist = cell(1,R0max+1);
for  qq = 1:R0max+vecdiff;
        
        try
        %load structure with trajectories
        if qq<(vecdiff+1)
            loadq = -1;
        else
            loadq = 0;

        end
        filenameload = strcat('data/stochsim_',country,'_R0',num2str(-vecdiff+loadq+qq),'_r1.mat');
        load(filenameload)
        
        %Set initial time of epidemic to zero for fitting
        ttimes = Infectdyncell.t(:,1)-min(Infectdyncell.t(:,1));
        
        %Estimate lambda for each trajectory
        numruns = size(Infectdyncell.t,2);  %this is 10000
        
        lambdadist{qq} = zeros(numruns,1);    %rate constant distribution
        zerotimedist{qq} = zeros(numruns,1);  %initial time of epidemic distribution
    
        %fit exponential to cumulative data
        for jj=1:numruns
            
            %fit line to log data
            [p,~] = fitdata_pois(ttimes,squeeze(sum(Infectdyncell.IRBDoft(:,:,jj),2)));
            
            %slope gives rate constant. save it, along with initial time
            lambdadist{qq}(jj,1) = p(1);
            zerotimedist{qq}(jj,1) = Infectdyncell.t(1,jj);
        end
    
    
    
    %Determine confidence intervals of rate constant distribution
    [data.lambda_low(qq),data.mlambda(qq),data.lambda_high(qq)] = diststats95(lambdadist{qq});
    
    %corresponding theoretical rate constant, useful for plotting lambda CI
    %plots

        end
end

data.xsupport = R0vec(1:R0max+vecdiff);
%Save full distribution
data.lambdadist = lambdadist;
data.zerotimedist = zerotimedist;

filenamesave = strcat('data/ci_',country,'_R0direct.mat');
save(filenamesave,'data')