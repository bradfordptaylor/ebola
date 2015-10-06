function plot_fig1
 
load('data/fig1_data')

rhodlist=.25;
%prune for dynamics for epidemics that reached the trigger and maintained
nonantime = data.starttime(~isnan(data.starttime));
nonaninfect = data.infectens(:,~isnan(data.starttime));

% 95% of data
sorttime = sort(nonantime);
hightime = sorttime(floor(.95*numel(nonantime)));
lowtime = sorttime(ceil(.05*numel(nonantime)));

%fit to exponential to get histogram
pens = zeros(1,size(nonaninfect,2));
for kk = 1:size(nonaninfect,2)
    currdyn = nonaninfect(:,kk);
    x = 0:length(currdyn)-1;
    [pens(kk),~] = fitdata_pois(x,currdyn);
end
sortpens = sort(pens);


%histogram figure

clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figures/fig1_hist';
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

[CIlow,CImed,CIhigh] = diststats95(1./pens)

%Shade region of CI in grey
hist(1./pens,linspace(min(1./pens),max(1./pens),50));
    
tmpt=text(max(1./pens)-.5*(max(1./pens)-min(1./pens)),300,sprintf('$\\tilde\\tau_c=%4.2f$ days',median(1./pens)));
  set(tmpt,'interpreter','latex','fontsize',20);
tmpt=text(max(1./pens)-.5*(max(1./pens)-min(1./pens)),200,sprintf('$\\rho_D=%4.2f$',rhodlist));
  set(tmpt,'interpreter','latex','fontsize',20);
xlabel('Characteristic time, $\hat{\tau}_c$ (days)','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Frequency','fontsize',20,'verticalalignment','bottom','interpreter','latex');
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

% automatic creation of postscript
psprintc(tmpfilename);

clear tmp*




%figure of dynamics across ensemble

%middle 95% of infected for every time point
dynlow = zeros(1,size(nonaninfect,1));
dynhigh = zeros(1,size(nonaninfect,1));
dynmed = zeros(1,size(nonaninfect,1));
lowidx = ceil(.05*numel(nonantime));
highidx = floor(.95*numel(nonantime));
for mm = 1:size(nonaninfect,1)
    dynatt = nonaninfect(mm,:);
    sortdyn = sort(dynatt);
    dynlow(mm) = sortdyn(lowidx);
    dynhigh(mm) = sortdyn(highidx);
    dynmed(mm) = median(dynatt);
end


xsupp = 0:42;   %x axis support
clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figures/fig1_dynwci';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

set(gcf,'DefaultLineMarkerSize',10);
set(gcf,'DefaultAxesLineWidth',2);
set(gcf,'PaperPositionMode','auto');



%Shade region of CI in grey
patch([xsupp fliplr(xsupp)],[dynlow fliplr(dynhigh)],[.8 .8 .8])

hold on
tmph = plot(xsupp,dynmed,'b-');
set(tmph,'linewidth',3,'MarkerSize',30);

tmph = plot(xsupp,dynhigh,'k-');
set(tmph,'linewidth',3,'MarkerSize',30);

tmph = plot(xsupp,dynlow,'k-');
set(tmph,'linewidth',3,'MarkerSize',30);

 
    hold off
    
%title(strcat(titlename),'fontsize',20,'interpreter','latex');
xlabel('Measurement time, days','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Cumulative cases','fontsize',20,'verticalalignment','bottom','interpreter','latex');
set(tmph,'linewidth',3);
set(gca,'fontsize',20);
xlim([0 42])

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





%part A -- Single Trajectories.
load('data/allcasecnts')

clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figures/fig1_traject';

tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);


set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');

set(gcf,'Position',[ 696          78         870        1237]);

maxt = 0;
for i=1:5,
  thetimes = transpose(0:(length(allcasecnts{i})-1));
  maxt = max([maxt length(thetimes)]);
end
for i=1:5,
  thetimes = transpose(0:(length(allcasecnts{i})-1));
  firstpt = thetimes(allcasecnts{i}>50);
  firstpt = min(firstpt);
  tmppos= [0.2 0.1+(0.15*(i-1)) 0.7 0.15];
  tmpa1 = axes('position',tmppos);
  Ioft = allcasecnts{i}(firstpt:end);
            
  %fit line to log data
  [p,~] = fitdata_pois(thetimes(firstpt:end),Ioft);
  tmph=plot(thetimes(1:firstpt),allcasecnts{i}(1:firstpt),'k.-');
  set(tmph,'linewidth',3,'color',[0.6 0.6 0.6]);
  hold on
  tmph=plot(thetimes(firstpt:end),allcasecnts{i}(firstpt:end),'k.-');
  set(tmph,'linewidth',3,'color','k');
  xlim([-10 25*ceil(maxt/25)]);
  set(gca,'xtick',[0:25:25*ceil(maxt/25)]);
  if (i~=1) 
    set(gca,'xticklabel',[]);
  else
    xlabel('Days, $t$','fontsize',24,'verticalalignment','top','interpreter','latex');
  end
  ylim([-50 550]);
  set(gca,'ytick',[0:100:500]);
  set(gca,'fontsize',20);
  tmpt=text(25,450,sprintf('$\\hat{\\tau}=%4.2f$ days',1/p(1)));
  set(tmpt,'interpreter','latex','fontsize',20);
  ylabel('Total Cases','fontsize',24,'verticalalignment','bottom','interpreter','latex');
end

tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');

psprintc(tmpfilenoname);

tmpt = pwd;
tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
text(1.05,.05,tmpnamememo,'Fontsize',6,'rotation',90);
datenamer(1.1,.05,90);

% automatic creation of postscript
psprintc(tmpfilename);

clear tmp*