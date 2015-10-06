function plot_errorhist

filenameload = strcat('data/ci_sim_cum_trig50.mat');
load(filenameload)
data_cum = data;

clear data

filenameload = strcat('data/ci_sim_inc_trig50.mat');
load(filenameload)
data_inc = data;

clear data


close all
%get histogram of standard errors of fits
thislamb = 7;
histpts = 200;
histlim_min = 0;
histlim_max = .01;
[yc,xc] = hist(data_cum.poisfitstderr{thislamb},linspace(histlim_min,histlim_max,histpts));
[yi,xi] = hist(data_inc.poisfitstderr{thislamb},linspace(histlim_min,histlim_max,histpts));

%plot
clf




eval(strcat('tmpfilename = ''figures/errorfitpois_large'';'))
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

%  tmph=plot(xc,yc/numel(data_cum.errordistlow{thislamb}), 'r--');
%  set(tmph,'linewidth',3);

plot(xc,yc/numel(data_cum.poisfitstderr{thislamb}),'LineWidth',3)
hold all
plot(xi,yi/numel(data_cum.poisfitstderr{thislamb}),'LineWidth',3)
legend('Cumulative','Incident')
legend boxoff
hold off
xlim([0 histlim_max])


xlabel('Fractional fitting error','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Probability Density','fontsize',20,'verticalalignment','bottom','interpreter','latex');
%title('Poisson Error','fontsize',20,'interpreter','latex');
set(gca,'fontsize',20);

tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');

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