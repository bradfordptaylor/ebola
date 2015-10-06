function [t,pop] = gill_onestep(t,pop,pars,numboxcar,rxneffect,kk)

probs = SEIRD(pop,pars,numboxcar,kk); %update rxn rates
alltaus = (1./probs).*log(1./rand(1,numboxcar+5));  %sample exponentially distributed times
[delt,rxnpos] = min(alltaus);   %minimum time chooses reaction
pop = pop + rxneffect(rxnpos,:);    %apply reaction to population
t = t+ delt;                        %update time