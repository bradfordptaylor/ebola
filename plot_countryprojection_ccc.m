function plot_countryprojection_ccc(country)


loadname = strcat('data/xsolve_',country);
load(loadname)

%get fits to data
switch country
    case 'guinea'
[p,s,cumdata,cumdata_all,xsupport] = fitcasedata_pois_data(country);
real_lambda = p(2);  % Exponential growth rate as fit to the data
titlename = 'Guinea';
    case 'liberia'
[p,s,cumdata,cumdata_all,xsupport] = fitcasedata_pois_data(country);
real_lambda = p(2);  % Exponential growth rate as fit to the data
        titlename = 'Liberia';
    case 'sleone'
[p,s,cumdata,cumdata_all,xsupport] = fitcasedata_pois_data(country);
real_lambda = p(2);  % Exponential growth rate as fit to the data
titlename = 'Sierra Leone';
end

%estimate tauc based on uniform distribution

%Standard simulation parameters
numboxcar = 2;
pars.Ti = 6;
pars.Te = 11;
pars.Td = 4;
pars.f = .7;
pars.ne = numboxcar;
%pars.rhod = rhodlist;
%pars.lambda = zeros(numruns,1);
%pars.R0 = R0;


delx = .005;
rhodvec = .1+delx:2*delx:.4-delx;
lambda1 = zeros(size(rhodvec));
lambda2 = zeros(size(rhodvec));
for kk=1:length(rhodvec)
    rhodtemp = rhodvec(kk);
    lambda1(kk) = R02lambda(xsolve(1),rhodtemp,pars);
    lambda2(kk) = R02lambda(xsolve(2),rhodtemp,pars);
end




numweeks = 8;
mlambda1 = mean(lambda1);
mlambda2 = mean(lambda2);
xtimes = 0:7:7*(length(cumdata)-1);
xtimesextend = 0:7:numweeks*7;
theoryfit = @(x) cumdata(end).*exp(p(2).*x);
theory1 = @(x) cumdata(end).*exp(mlambda1.*x);
theory2 = @(x) cumdata(end).*exp(mlambda2.*x);



clf;
% automatically create postscript whenever
% figure is drawn
eval(strcat('tmpfilename = ''figures/figcountry_dataproj_',country,''';'))
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');



%Shade region of CI in grey
patch([xtimes(end)+xtimesextend fliplr(xtimes(end)+xtimesextend)],[theory1(xtimesextend) fliplr(theory2(xtimesextend))],[.8 .8 .8])

hold on

tmph = plot(xtimes,cumdata,'k.');
set(tmph,'linewidth',3,'MarkerSize',30);

tmph = plot(xtimes(end)+xtimesextend(2:end),cumdata_all(numel(xtimes)+1:numel(xtimes)+numel(xtimesextend)-1),'gd','MarkerFaceColor','g');
set(tmph,'linewidth',3,'MarkerSize',10);

tmph = plot(xtimes(end)+xtimesextend,theoryfit(xtimesextend),'b-');
set(tmph,'linewidth',3,'MarkerSize',30);

tmph = plot(xtimes(end)+xtimesextend,theory1(xtimesextend),'k-');
set(tmph,'linewidth',3,'MarkerSize',30);

tmph = plot(xtimes(end)+xtimesextend,theory2(xtimesextend),'k-');
set(tmph,'linewidth',3,'MarkerSize',30);

xlim([0 max(xtimes(end)+xtimesextend)])
%  set(tmph,'linewidth',3);
%    tmph=plot(xlim, [1./real_lambda 1./real_lambda], 'r-');
%    set(tmph,'linewidth',3);
    hold off
%title(strcat(titlename),'fontsize',20,'interpreter','latex');
xlabel('Time (days)','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Cumulative case count','fontsize',20,'verticalalignment','bottom','interpreter','latex');
title(titlename,'fontsize',20,'interpreter','latex');
set(tmph,'linewidth',3);
set(gca,'fontsize',20,'yscale','linear');

tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');

% first two points are normalized x, y positions
% text(,,'','Fontsize',14);

% automatic creation of postscript
% without name/date
psprintc(tmpfilenoname);
psprint(tmpfilebwname);

tmpt = pwd;
tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
text(1.05,.05,tmpnamememo,'Fontsize',6,'rotation',90);
datenamer(1.1,.05,90);
% datename(.5,.05);
% datename2(.5,.05); % 2 rows

% automatic creation of postscript
psprintc(tmpfilename);

% set following on if zooming of 
% plots is required
% may need to get legend up as well
%axes(tmpa1)
%axes(tmplh)
clear tmp*