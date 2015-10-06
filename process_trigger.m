%%
%load data to get parameters of data
numsamps = 100;
filepath = strcat('~/Desktop/ebola_data_final/sims/simdata_final_L',num2str(1),'_r',num2str(1),'.mat');
load(filepath)
runsinthing = size(Infectdyncell.popdyn,3);
%preallocate
sampdays = 42;
alltimes = size(Infectdyncell.popdyn,1);
allIRBD = cell(1,20);
lambdalist = zeros(1,20);

Iveclist = 10:10:50;
for ll=1:20
    for mm = 1:length(Iveclist);
        allIRBD{ll}{mm} = zeros(sampdays+1+1,4,runsinthing*numsamps);
        count = 0;
    end
    
    for rr = 1:numsamps
            filepath = strcat('~/Desktop/ebola_data_final/sims/simdata_final_L',num2str(ll),'_r',num2str(rr),'.mat');
            load(filepath)
            count = count+1;
            for pp = 1:size(Infectdyncell.popdyn,3)
            for nn=1:length(Iveclist)
                currcum = sum(Infectdyncell.popdyn(:,3:end,pp),2);
                clipdyn = find(currcum>=Iveclist(nn),1);
                allIRBD{ll}{nn}(:,:,(count-1)*runsinthing+pp) = Infectdyncell.popdyn(clipdyn-1:clipdyn+sampdays,3:end,pp);   %all dynamics of IRBD classes
            end
            end

            
            
    end
   lambdalist(ll) = Infectdyncell.pars.lambda; %to use as x-support in later plots
end
save('data/sim_lambdalist','lambdalist')


for mm = 1:length(Iveclist)
Imin = Iveclist(mm);
cumdyntrigger = cell(1,20);
incdyntrigger = cell(1,20);

cumdyn = cell(1,20);
incdyn = cell(1,20);
for ll=1:20
    cumdyn{ll} = squeeze(sum(allIRBD{ll}{mm},2));
    incdyn{ll} = zeros(size(cumdyn{ll}));
    incdyn{ll}(1,:) = cumdyn{ll}(1,:);
    incdyn{ll}(2:end,:) = diff(cumdyn{ll});
end


for ll=1:20
    cumdyntrigger{ll} = zeros(43,size(cumdyn{ll},2));
    incdyntrigger{ll} = zeros(43,size(cumdyn{ll},2));
    for kk = 1:size(cumdyn{ll},2)
        lastbitcum = cumdyn{ll}(cumdyn{ll}(:,kk)>=Imin,kk);
        lastbitinc = incdyn{ll}(cumdyn{ll}(:,kk)>=Imin,kk);
    	cumdyntrigger{ll}(:,kk) = lastbitcum(1:43);
        incdyntrigger{ll}(:,kk) = lastbitinc(1:43);
    end
end
filesavecum = strcat('data/sim_cum_trig',num2str(Imin));
filesaveinc = strcat('data/sim_inc_trig',num2str(Imin));
save(filesavecum,'cumdyntrigger')
save(filesaveinc,'incdyntrigger')
end