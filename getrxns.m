function rxneffect = getrxns(numboxcar)

%Effect on population vector for each process
rxneffect = zeros(numboxcar+5,numboxcar+5); %initialize
rxneffect(1,:) = [-1 1 zeros(1,numboxcar-1) 0 0 0 0];   %Infect from infectious
rxneffect(2,:) = [-1 1 zeros(1,numboxcar-1) 0 0 0 0];   %Infect from deceased

for jj = 1:numboxcar
    rxneffect(jj+2,:) = [zeros(1,jj) -1 1 zeros(1,numboxcar-jj) 0 0 0]; %Progress through exposed class
end
rxneffect(numboxcar+3,:) = [0 zeros(1,numboxcar) -1 1 0 0]; %Infectious recovers
rxneffect(numboxcar+4,:) = [0 zeros(1,numboxcar) -1 0 1 0]; %Infectious dies
rxneffect(numboxcar+5,:) = [0 zeros(1,numboxcar) 0 0 -1 1]; %Dead buried