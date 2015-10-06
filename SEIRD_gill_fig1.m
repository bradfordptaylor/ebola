function [infectdyn,iscrash,starttime] = SEIRD_gill_fig1(x0, Imin,pars,numboxcar)

if nargin==3 || numboxcar<1
    numboxcar=1;
end

tf = 42;    %total time observed after trigger
tsamp_rel = 0:tf;
tsteps = length(tsamp_rel);
tsamp_rel = [tsamp_rel 1000*tf];   %pad end with a large number for while loop

rxneffect = getrxns(numboxcar);

%Start popution with single infectious individual
pop = [x0-1 zeros(1,numboxcar) 1 0 0 0];   %[S Evec I R D B]
totclass = numboxcar + 5;   %number of state variables

%gillespie
t=0;    %set time for beginning of epidemic
%priorI = 1;
priorcum = 1;   %initially 1 cumulative infection
nextsamp = 1;
%run dynamics until trigger number of cumulative infected, Imin, is reached
while  priorcum<Imin && sum(pop(2:end-3))+pop(end-1)>0
    prevdata = pop;
    [t,pop] = gill_onestep(t,pop,pars,numboxcar,rxneffect);
    if t>=nextsamp
        priorcum = (prevdata(1,totclass-3)+prevdata(1,totclass-2)+prevdata(1,totclass-1)+prevdata(1,totclass));
        nextsamp = ceil(t);
    end
end
%label epidemics that dissipate with NaNs else prep to sample case data
if priorcum<Imin
    iscrash = 1;
    starttime = NaN;
    infectdyn = NaN*ones(tsteps,1);
else
    starttime = floor(t);
    t0actual = t;
    popdyn = zeros(length(tsamp_rel)-1,totclass);
    t = t0actual-floor(t0actual);
    count = 1;
    %set first data point of case data
    popdyn(count,:) = prevdata;
    count = count+1;
    %sample case data
    while t<tf && sum(pop(2:end-3))+pop(end-1)>0
        prevdata = pop;
        [t,pop] = gill_onestep(t,pop,pars,numboxcar,rxneffect);
        %sample next integer time point
        while t>=tsamp_rel(count)
            popdyn(count,:) = prevdata;
            count=count+1;
        end
        
    end

    if t<tf %&& sum(pop(2:end-3))+pop(end-1)==0
        iscrash = 1;
        starttime = NaN;
        infectdyn = NaN*ones(tsteps,1);
    else
        iscrash=0;
        infectdyn = sum(popdyn(:,numboxcar+2:numboxcar+5),2);   %cumulative infection dynamics

    end
    
    
end
end
