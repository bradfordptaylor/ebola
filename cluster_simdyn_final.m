function cluster_simdyn_final(lambdan,runseq)

numruns = 100;

%number of exposed classes in simulation
numboxcar = 1;
totclass = numboxcar + 5;

%define params, as in figure 3
pars.Ti = 6;
pars.Te = 11;
pars.Td = 4;
pars.f = .7;
lambdalist = 1./linspace(14,33,20);%[1/14 1/21 1/28];



N = 10^6; %inital population
Imax = 50;   %last trigger size to consider
rhod0 = .25;
pars.rhod = rhod0;

lambda0 = lambdalist(lambdan);
pars.lambda = lambda0;

tf = 500;   %set a long time to catch slow epidemics
tsamp = transpose(0:tf);
        
%PreAllocate
allt = zeros(tf+1,numruns);
allpopdyn = zeros(tf+1,totclass,numruns);

%Derive other simulation parameters 
if rhod0 ==0
	pars.betad = 0;
	pars.betai = ((1+pars.Ti*lambda0)*(1+(pars.Te*lambda0/numboxcar))^numboxcar)/pars.Ti;
else
	pars.betad = -((1 + pars.Td*lambda0)*(1 + pars.Ti*lambda0)*rhod0*(1 + (pars.Te*lambda0)/numboxcar)^numboxcar)...
        /(pars.f*pars.Td*(-1 + pars.Td*lambda0*(-1 + rhod0)));
pars.betai = -((1 + pars.Td*lambda0)*(1 + pars.Ti*lambda0)*(1-rhod0)*(1 + (pars.Te*lambda0)/numboxcar)^numboxcar)...
    /(pars.Ti*(-1+pars.Td*lambda0*(rhod0-1)));
end
        
        
%loop over runs to obtain trajectories of cumulative infected
parfor jj=1:numruns

    %Simulate
    tempdata = SEIRD_gill_fulldyn(N,Imax,tsamp,pars,numboxcar);
	allpopdyn(:,:,jj) = tempdata.popdyn;
    %Record time points in simulation corresponding to census times
	allt(:,jj) = tempdata.t;
    %alltfinal(1,jj) = tempdata.tfinal;
end

%Place data in one structure, Infectdyncell
Infectdyncell.popdyn = allpopdyn;
Infectdyncell.t = allt;
Infectdyncell.pars = pars;
%Infectdyncell.tfinal = alltfinal;
        

%save as data/stochsim_guinea_rhod1_L1_r1 for first rhod, lambda, and
%runcase. Iterate through and save.
savename = strcat('data/simdata_final'...
    ,'_L',num2str(lambdan),'_r',num2str(runseq));
save(savename,'Infectdyncell')
