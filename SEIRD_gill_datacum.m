function data = SEIRD_gill_datacum(x0, Imin,tsamp,pars,numboxcar,jj)
 
startday = min(tsamp);
tsamp_rel = tsamp-startday; %set times to start at 0
tf = max(tsamp_rel);
tsamp_rel = [tsamp_rel; 1000*tf]; %pad with large number for while loop

rxneffect = getrxns(numboxcar);

%initial population is 1 infected
pop = [x0-1 zeros(1,numboxcar) 1 0 0 0];   %[S Evec I R D B]
totclass = numboxcar + 5;

%run until trigger
t=0;
priorI = 1;
priorcum = 1;
precount = 1;
while  priorcum<Imin && sum(pop(2:end-3))+pop(end-1)>0
    prevdata = pop;
    [t,pop] = gill_onestep(t,pop,pars,numboxcar,rxneffect,jj);
    while t>=(precount*7)
        priorI = sum(prevdata(end-3:end))-priorcum;
        priorcum = sum(prevdata(end-3:end));
        precount = precount+1;
    end
end
t0actual = t;



if sum(pop(2:end-3))+pop(1,end-1)==0
    %If no epidemic, rerun
    data = SEIRD_gill_datacum(x0, Imin,tsamp,pars,numboxcar,jj);
else
%else continue


popdyn = zeros(length(tsamp),totclass);
%t=0;
t = mod(t0actual,7);    %weekly
popdyn(1,:) = prevdata;
count = 2;
%sample case data
while t<tf && sum(pop(2:end-3))+pop(end-1)>0
    prevdata = pop;
    [t,pop] = gill_onestep(t,pop,pars,numboxcar,rxneffect,jj);

    while t>=tsamp_rel(count)
        popdyn(count,:) = prevdata;
        count=count+1;
    end


    
end

 if t<tf %&& sum(pop(2:end-3))+pop(end-1)==0
     %If no epidemic, rerun
     data = SEIRD_gill_datacum(x0,Imin,tsamp,pars,numboxcar,jj);
 else
     %export data

data.t = tsamp_rel(1:end-1);
data.Ioft = popdyn(:,numboxcar+2);
data.Roft = popdyn(:,numboxcar+3);
data.Boft = popdyn(:,numboxcar+4);
data.Doft = popdyn(:,numboxcar+5);
data.priorI = priorI;

 end


 end
end
