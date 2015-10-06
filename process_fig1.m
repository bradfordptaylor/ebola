function process_fig1

%Parameters for Sims

numtraj = 10^4;	%number of trajectories in ensemble
lambda0 = 1/21; %Theoretical characteristic time of growth
numboxcar = 1;  %number of exposed classes in simulation
Imin=50;        %Trigger population
N = 10^6;       %Total population


%Model parameters
rhod0 = .25; %fraction of transmissions from deceased for fig 1
pars.Ti = 6;
pars.Te = 11;
pars.Td = 4;
pars.f = .7;
pars.ne = numboxcar;

pars.rhod = rhod0;

pars.lambda = lambda0;
        


%Derive other simulation parameters 
if rhod0 ==0
    pars.betad = 0;
    pars.betai = ((1+pars.Ti*lambda0)*(1+(pars.Te*lambda0/numboxcar))^numboxcar)/pars.Ti;
else
    pars.betad = -((1 + pars.Td*lambda0)*(1 + pars.Ti*lambda0)*rhod0*(1 +...
        (pars.Te*lambda0)/numboxcar)^numboxcar)/(pars.f*pars.Td*(-1 +...
        pars.Td*lambda0*(-1 + rhod0)));
    pars.betai = -((1 + pars.Td*lambda0)*(1 + pars.Ti*lambda0)*(1-rhod0)*...
        (1 + (pars.Te*lambda0)/numboxcar)^numboxcar)...
        /(pars.Ti*(-1+pars.Td*lambda0*(rhod0-1)));
end
        
%loop for stats
%loop over runs to obtain trajectories of cumulative infected
parfor jj=1:numtraj
    [infectens(:,jj),iscrash(jj),starttime(jj)] = SEIRD_gill_fig1(N,Imin,pars,numboxcar);
end

data.infectens = infectens;
data.iscrash = iscrash;
data.starttime = starttime;

save('data/fig1_data','data')


%single trajectories for part A      
allcasecnts = cell(1,5);

for jj=1:5
infectens = SEIRD_gill_fig1_traj(N,Imin,pars,numboxcar);
thetimes = transpose(0:(length(infectens)-1));
firstpt = min(thetimes(infectens>50));
allcasecnts{jj} = infectens(1:firstpt+43);
end
save('data/allcasecnts','allcasecnts')