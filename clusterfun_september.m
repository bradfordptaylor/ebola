function clusterfun_september(country,R0seq,runseq,numruns)
%
% Case fitting

rng(runseq)
R0list = 1.1:.1:2.5;
try
R0 = R0list(R0seq);
catch
    %extend
    R0list2 = [1.05 1.025 1.0125];
    R0 = R0list2(abs(R0seq));
end
%sample from prior for rho_D
rhodlist = .3.*rand(numruns,1)+.1;

%number of exposed classes in simulation
numboxcar = 2;

%switch between countries to define the
%country-specific parameters for the simulations, case 1 (guinea), case 2 is liberia,
%and case 3 is sierra leone


%census times matching the case data
switch country
    case 1
        tempcount = 'guinea';
        [~,~,data_cum,~,~,~,~] = fitcasedata_pois_data(tempcount);
        data_cum= data_cum(data_cum>=50);
        Imin = data_cum(1);
        %census times matching the case data
        tsamp = transpose(0:7:7*(length(data_cum)-1));


        %Total population (guinea specific here)
        N = 11.8*10^6;
        
    case 2
        
        tempcount = 'liberia';
        [~,~,data_cum,~,~,~,~] = fitcasedata_pois_data(tempcount);
        data_cum= data_cum(data_cum>=50);
        Imin = data_cum(1);
        %census times matching the case data
        tsamp = transpose(0:7:7*(length(data_cum)-1));

        %Total population (guinea specific here)
        N =4.3*10^6;
        
    case 3
        
        tempcount = 'sleone';
        [~,~,data_cum,~,~,~,~] = fitcasedata_pois_data(tempcount);
        data_cum= data_cum(data_cum>=50);
        Imin = data_cum(1);
        %census times matching the case data
        tsamp = transpose(0:7:7*(length(data_cum)-1));

        %Total population (guinea specific here)
        N=6.1*10^6;
end


%Standard simulation parameters
pars.Ti = 6;
pars.Te = 11;
pars.Td = 4;
pars.f = .7;
pars.ne = numboxcar;
pars.rhod = rhodlist;
pars.lambda = zeros(numruns,1);
pars.R0 = R0;

for kk=1:numruns
    rhod = rhodlist(kk);
    lambda0 = R02lambda(R0,rhod,pars);
    pars.lambda(kk) = lambda0;
end
pars.betad = -((1 + pars.Td.*pars.lambda).*(1 + pars.Ti.*pars.lambda).*pars.rhod.*(1 + (pars.Te.*pars.lambda)./numboxcar).^numboxcar)...
                        ./(pars.f.*pars.Td.*(-1 + pars.Td.*pars.lambda.*(-1 + pars.rhod)));
pars.betai = -((1 + pars.Td.*pars.lambda).*(1 + pars.Ti.*pars.lambda).*(1-pars.rhod).*(1 + (pars.Te.*pars.lambda)./numboxcar).^numboxcar)...
                        ./(pars.Ti.*(-1+pars.Td.*pars.lambda.*(pars.rhod-1)));
        
allIRBDoft = zeros(length(tsamp),4,numruns);
allt = zeros(length(tsamp),numruns);
allpriorI = zeros(1,numruns);
        
        
        %loop over runs to obtain trajectories of cumulative infected
        parfor jj=1:numruns

                    
                %Simulate
                tempdata = SEIRD_gill_datacum(N,Imin,tsamp,pars,numboxcar,jj);
                %If the last element is 0 then the epidemic ended before
                %the final census time, throw out this result and rerun
            %Record Infected, Recovered, Buried, and dead. Their sum is
            %cumulative infected
            allIRBDoft(:,:,jj) = [tempdata.Ioft tempdata.Roft tempdata.Boft tempdata.Doft];
            %Record time points in simulation corresponding to census times
            allt(:,jj) = tempdata.t;
            allpriorI(jj) = tempdata.priorI;
        end
        
        %Place data in one structure, Infectdyncell
        data.IRBDoft = allIRBDoft;
        data.t = allt;
        data.pars = pars;
        data.rho = rhodlist;
        data.priorI = allpriorI;
        Infectdyncell = data;
        


%save data
savename = strcat('data/stochsim_',tempcount,'_R0',num2str(R0seq)...
    ,'_r',num2str(runseq));
save(savename,'Infectdyncell')



end
