function infectdyn = SEIRD_gill_fig1_traj(x0, Imin,pars,numboxcar)

totalt = 225;
tsteps = 42;

if nargin==3 || numboxcar<1
    numboxcar=1;
end
tsamp_tot =0:totalt;
tsamp_tot =[tsamp_tot 10*totalt];

rxneffect = getrxns(numboxcar);

%x0 = 1
pop = [x0-1 zeros(1,numboxcar) 1 0 0 0];   %[S Evec I R D B]
totclass = numboxcar + 5;

%gillespie

popdyn = zeros(length(tsamp_tot)-1,totclass);

t = 0;
count = 1;

%sample case data
while t<totalt && sum(pop(2:end-3))+pop(end-1)>0
    prevdata = pop;
    [t,pop] = gill_onestep(t,pop,pars,numboxcar,rxneffect);
    while t>=tsamp_tot(count)
        popdyn(count,:) = prevdata;
        count=count+1;
    end
end

%if epidemic dissipates rerun
 if t<totalt 
    infectdyn = SEIRD_gill_fig1_traj(x0, Imin,pars,numboxcar);
 else
     infectdyn = sum(popdyn(:,numboxcar+2:numboxcar+5),2);  %cumulative infected

 end


end
