function data = SEIRD_gill_fulldyn(x0, Imax,tsamp,pars,numboxcar)

if nargin==4 || numboxcar<1
    numboxcar=1;
end
tf = max(tsamp);

rxneffect = getrxns(numboxcar);

%population at t=0
pop = [x0-1 zeros(1,numboxcar) 1 0 0 0];   %[S Evec I R D B]
totclass = numboxcar + 5;

%gillespie

t=0;

popdyn = zeros(length(tsamp),totclass);
%t=0;
count = 1;
cumI = 1;
while t<tf && sum(pop(2:end-3))+pop(end-1)>0 && cumI<Imax
    prevdata = pop;
    [t,pop] = gill_onestep(t,pop,pars,numboxcar,rxneffect,1);    
    while t>=tsamp(count)
        popdyn(count,:) = prevdata;
        count=count+1;
        cumI = sum(prevdata(end-3:end));
    end

end
%get out of loop if t>tf, epidemic ends, or cumulative is bigger than Imax

 if sum(pop(2:end-3))+pop(end-1)==0 %Epidemic ends prematurely
     %If no epidemic, rerun
     data = SEIRD_gill_fulldyn(x0,Imax,tsamp,pars,numboxcar);
 else
     %run for 42 more points
     currfill = 0;
     
     while t<tf && sum(pop(2:end-3))+pop(end-1)>0 && currfill<42
    prevdata = pop;
    [t,pop] = gill_onestep(t,pop,pars,numboxcar,rxneffect,1);    
    while t>=tsamp(count)
        popdyn(count,:) = prevdata;
        count=count+1;
        currfill = currfill+1;
    end

     end
     %if the epidemic ended or t>tf rerun
     if currfill~=42
         data = SEIRD_gill_fulldyn(x0,Imax,tsamp,pars,numboxcar);
     else

          %data.t = t0actual+tsamp_rel;
data.t = tsamp;
data.popdyn = popdyn;
%data.tfinal = t-delt;

     end
     %data = SEIRD_gill(x0, Imin, tf,numpatch,pars,numboxcar); old

 end
