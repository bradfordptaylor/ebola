function [r,stderr,data_cum,data_cumall,data_inc,data_incall,xsupport] = fitcasedata_pois_data(country)

%filepath = strcat('data/',country,'_121714.csv');
filepath = 'data/casedata_countries_gls.csv';
data_incident = csvread(filepath,1);
data_early = data_incident(1:36,:);
data_all = data_incident;
xsupport = data_all(:,2);
switch country
    case 'guinea'
data_cum = cumsum(data_early(:,3));
data_cumall = cumsum(data_all(:,3));
data_inc = data_early(:,3);
data_incall = data_all(:,3);
datacumlog = data_cumall>=50;
data_cumall = data_cumall(datacumlog);
data_incall = data_incall(datacumlog);
xsupport = xsupport(datacumlog);
data_inc= data_inc(data_cum>=50);
data_cum= data_cum(data_cum>=50);

[r,stderr] = fitdata_pois_data(0:7:7*(length(data_cum)-1),data_cum');
    case 'liberia'
data_cum = cumsum(data_early(:,4));
data_cumall = cumsum(data_all(:,4));
data_inc = data_early(:,4);
data_incall = data_all(:,4);
datacumlog = data_cumall>=50;
data_cumall = data_cumall(datacumlog);
data_incall = data_incall(datacumlog);
xsupport = xsupport(datacumlog);
data_inc= data_inc(data_cum>=50);
        data_cum= data_cum(data_cum>=50);
        
[r,stderr] = fitdata_pois_data(0:7:7*(length(data_cum)-1),data_cum');
    case 'sleone'
data_cum = cumsum(data_early(:,5));
data_cumall = cumsum(data_all(:,5));
data_inc = data_early(:,5);
data_incall = data_all(:,5);
datacumlog = data_cumall>=50;
data_cumall = data_cumall(datacumlog);
data_incall = data_incall(datacumlog);
xsupport = xsupport(datacumlog);
data_inc= data_inc(data_cum>=50);
        data_cum= data_cum(data_cum>=50);
[r,stderr] = fitdata_pois_data(0:7:7*(length(data_cum)-1),data_cum');
end