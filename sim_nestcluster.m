function sim_nestcluster(thislambda,startidx)

%startvec = [1 25 50 76];
startvec = 1:10:100;
startnum = startvec(startidx);

%startvec = [46 43 42 59 56 53 42 41 48 48 41 47 44 45 42 37 39 38 38 37];
%startvec = startnum*ones(size(startvec));
totalruns = 100;
for mm=startnum:totalruns
    
    try
        load(strcat('data/simdata_final_L',num2str(thislambda),'_r',num2str(mm),'.mat'));
        clear Infectdyncell
    catch
    rng(thislambda+mm*21)
    cluster_simdyn_final(thislambda,mm)
    end
end