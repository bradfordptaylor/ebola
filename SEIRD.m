function prob = SEIRD(pop,pars,numboxcar,kk)

%define Reactions based on the following dynamical system:

%Sdot = -betai*S.*I./N-betad*S.*D./N
%Edot = betai*S.*I./N +betad*S.*D./N - E./Te
%Idot = E/Te - I/Ti
%Rdot = (1-f)*I./Ti
%Ddot = f*I./Ti - D/Td
%Bdot = D/Td

totclass = 5+numboxcar;
prob = zeros(1,numboxcar+5);    %initialize

prob(1) = pars.betai(kk)*pop(1)*pop(2+numboxcar)./sum(pop(1:end-2));    %Infect from infectious


prob(2) = pars.betad(kk)*pop(1)*pop(totclass-1)./sum(pop(1:end-2));     %Infect from deceased


for jj = 1:numboxcar
    prob(jj+2) = (numboxcar/pars.Te)*pop(jj+1);     %Progress through exposed class
end

prob(numboxcar+3) = (1/pars.Ti)*(1-pars.f)*pop(totclass-3); %Infectious recovers


prob(numboxcar+4) = (1/pars.Ti)*(pars.f)*pop(totclass-3);   %Infectious dies


prob(numboxcar+5) = (1/pars.Td)*pop(totclass-1);    %Dead buried
