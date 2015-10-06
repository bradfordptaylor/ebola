function plot_sim_distcomp

%load data
filenameload = strcat('data/ci_sim_cum_trig50.mat');
load(filenameload)
data_cum = data;

clear data

filenameload = strcat('data/ci_sim_inc_trig50.mat');
load(filenameload)
data_inc = data;

clear data

%look at hist for a given lambda
lambdachoice = 7; %corresponds to tau_c = 20
histpts = 50;

histlim_min = min([min(1./data_cum.poislambdadist{lambdachoice}) min(1./data_inc.poislambdadist{lambdachoice})]);
histlim_max = max([max(1./data_cum.poislambdadist{lambdachoice}) max(1./data_inc.poislambdadist{lambdachoice})]);
histlim_max = 40;   %set because of long tail
cleancum = data_cum.poislambdadist{lambdachoice}((1./data_cum.poislambdadist{lambdachoice})<=histlim_max);
cleaninc = data_inc.poislambdadist{lambdachoice}((1./data_inc.poislambdadist{lambdachoice})<=histlim_max);

%plot

clf
eval(strcat('tmpfilename = ''figures/fig_cuminc_poisdistcompare'';'))
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

[ys,xs] = hist(1./cleancum,linspace(histlim_min,histlim_max,histpts));
ys = ys/numel(1./cleancum);
h1 = plot(xs,ys,'LineWidth',3)
%alpha(.5)
 hold all
 %patch([xs fliplr(xs)],[zeros(size(ys)) fliplr(ys)],[.2 .8 .8])
 %alpha(.5)
 [ys2,xs2] = hist(1./cleaninc,linspace(histlim_min,histlim_max,histpts));
 ys2 = ys2/numel(1./cleaninc);
 h2 = plot(xs2,ys2,'LineWidth',3)
  %patch([xs2 fliplr(xs2)],[zeros(size(ys2)) fliplr(ys2)],[.8 .5 .2])
 %alpha(.5)
 plot([20 20],[0 1.1*max([ys ys2])],'k--','LineWidth',3)
% alpha(.5)
 hold off
 legend([h1 h2],{'cumulative','incident'})
 legend boxoff
 xlim([min([xs xs2]) max([xs xs2])])
 ylim([0 1.1*max([ys ys2])])
 xlabel('Measured $\tau_c$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Probability density','fontsize',20,'verticalalignment','bottom','interpreter','latex');
title(strcat('Theoretical $\tau_c$ =',num2str(20)),'fontsize',20,'interpreter','latex');
set(gca,'fontsize',20);

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