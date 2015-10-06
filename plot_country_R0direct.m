function plot_country_R0direct(country)

%load from process file
filenameload = strcat('data/ci_',country,'_R0direct.mat');
load(filenameload)



%measure growth rate from case data
[r,~,~,~,~,~,~] = fitcasedata_pois_data(country);
real_lambda = r(2);  % Exponential growth rate as fit to the data
switch country
    case 'guinea'


titlename = 'Guinea';
    case 'liberia'

        titlename = 'Liberia';
    case 'sleone'
titlename = 'Sierra Leone';

end

%find intersection of measured with simulated CI

    aboveptlow = find(diff(sign(1./data.lambda_low-1./real_lambda)));
    abovepthigh = find(diff(sign(1./data.lambda_high-1./real_lambda)));
    yone = [data.lambda_low(aboveptlow) data.lambda_high(abovepthigh)].^(-1);
    ytwo = [data.lambda_low(aboveptlow+1) data.lambda_high(abovepthigh+1)].^(-1);
    xone = [data.xsupport(aboveptlow) data.xsupport(abovepthigh)];
    xtwo = [data.xsupport(aboveptlow+1) data.xsupport(abovepthigh+1)];
    xsolve = ((1./real_lambda-yone).*(xtwo-xone)./(ytwo-yone))+xone;

    savename = strcat('data/xsolve_',country);
    save(savename,'xsolve')
 
 
 
 
 

%plot
clf;
% automatically create postscript whenever
% figure is drawn
eval(strcat('tmpfilename = ''figures/figcountry_ciR0_',country,'_large'';'))
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

data.xsupport = data.xsupport(data.lambda_low~=0);
data.lambda_low = data.lambda_low(data.lambda_low~=0);
data.lambda_high = data.lambda_high(data.lambda_high~=0);
data.mlambda = data.mlambda(data.mlambda~=0);

%Shade region of CI in grey
patch([data.xsupport fliplr(data.xsupport)],[1./data.lambda_low fliplr(1./data.lambda_high)],[.8 .8 .8])

hold on
tmph = plot(data.xsupport,1./data.mlambda,'b.-');
set(tmph,'linewidth',3,'MarkerSize',30);

 tmph=plot(data.xsupport,1./data.lambda_low,'k.-');
 set(tmph,'linewidth',3,'MarkerSize',30);
 tmph=plot(data.xsupport,1./data.lambda_high,'k.-');
 set(tmph,'linewidth',3,'MarkerSize',30);

% ROmed = R0(real_lambda)
R0CI = xsolve
 %if ~strcmp(country,'guinea')
 tmph=plot(xsolve, [1./real_lambda 1./real_lambda], 'r-');
  set(tmph,'linewidth',3);
  tmph=plot([xsolve(1) xsolve(1)], [0 1./real_lambda], 'r--');
  set(tmph,'linewidth',3);
  tmph=plot([xsolve(2) xsolve(2)], [0 1./real_lambda], 'r--');
  set(tmph,'linewidth',3);
  %else
  %    tmph=plot([min(data.xsupport) xsolve], [1./real_lambda 1./real_lambda], 'r-');
  % set(tmph,'linewidth',3);
  %      tmph=plot([xsolve xsolve], [0 1./real_lambda], 'r--');
  % set(tmph,'linewidth',3); 
  %end


%Plot style

if numel(xsolve)==2
xlim([floor(10*min(xsolve))./10-.05 ceil(10*max(xsolve))./10+.05])
else
    ylim([min(1./data.lambda_high) 1.01*max(1./data.lambda_low)])
end
%  set(tmph,'linewidth',3);
%    tmph=plot(xlim, [1./real_lambda 1./real_lambda], 'r-');
%    set(tmph,'linewidth',3);
    hold off
title(strcat(titlename),'fontsize',20,'interpreter','latex');
xlabel('$R_0$ theory','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('$\tau_c$ measured','fontsize',20,'verticalalignment','bottom','interpreter','latex');
set(tmph,'linewidth',3);
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