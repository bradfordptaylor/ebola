%%
clf
clear all
close all
%plot sims using incident data- R0 CI
rhod = .25;


%Parameters for R0 in terms of lambda via generating function approache
%following Wallinga and Lipstitch

n_E    = 1;        % Number of stages for Gamma
info.T_E    = 11;       % Days exposed, on average
b_E    = n_E/info.T_E;  % Control parameter for Gamma
info.T_I    = 6;        % Number of days infectious
gamma    = 1/info.T_I;% Control parameter for I
info.T_D    = 4;        % Number of days until burial
chi    = 1./info.T_D;% Control parameter for D


ME = @(lambda) (b_E./(b_E-lambda)).^n_E;
MI = @(lambda) gamma./(gamma-lambda);
MD = @(lambda) chi./(chi-lambda);
   M = @(lambda) (1-rhod).*ME(lambda).*MI(lambda)+rhod.*ME(lambda).*MI(lambda).*MD(lambda);
R0 = @(lambda) 1./(M(-lambda));



%load CI boundary data from previous plot
filenameload = strcat('data/ci_sim_inc_trig50.mat');
load(filenameload)
loadname = strcat('data/xsolve_simci_inc');
load(loadname)
loadname = strcat('data/lambsolve_simci_inc');
load(loadname)
loadname = strcat('data/interpvec_simci_inc');
load(loadname)




real_lambda = 1/20;  % mock data





%Find intersection with R0 CI
    
     pt1 = xsolve(1);   %lambda_1 boundary
     pt2 = xsolve(2);   %lambda_2 boundary
    
    %Upper boundary intersection
    abovepthigh = find(diff(sign(lambsolve(:,1)-pt1)));
    yone = R0(lambsolve(abovepthigh,1).^-1);
    ytwo = R0(lambsolve(abovepthigh+1,1).^-1);
    xone = interpvec(1,abovepthigh);
    xtwo = interpvec(1,abovepthigh+1);
    ysolve = (ytwo-yone).*(pt1-xone)./(xtwo-xone)+yone;
        
    %Lower boundary intersection
    aboveptlow = find(diff(sign(lambsolve(:,2)-pt2)));
    yone = R0(lambsolve(aboveptlow,2).^-1);
    ytwo = R0(lambsolve(aboveptlow+1,2).^-1);
    xone = interpvec(2,aboveptlow);
    xtwo = interpvec(2,aboveptlow+1);
    ysolve2 = (ytwo-yone).*(pt2-xone)./(xtwo-xone)+yone;

clf;
% automatically create postscript whenever
% figure is drawn


eval(strcat('tmpfilename = ''figures/fig_simciR0_inc_large'';'))
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





%Shade CI range in grey
patch([interpvec(1,:) fliplr(interpvec(2,:))],R0([lambsolve(:,1)'.^-1 fliplr(lambsolve(:,2)'.^-1)]),[.8 .8 .8]);%R0([data{kk}.lambda_low fliplr(data{kk}.lambda_high)]),[.8 .8 .8])
      
hold on


tmph=plot(interpvec(3,:),R0(lambsolve(:,3)'.^-1),'b-');    %median
set(tmph,'linewidth',3,'MarkerSize',30);

tmph=plot(interpvec(1,:),R0(lambsolve(:,1)'.^-1),'k-');    %top CI
set(tmph,'linewidth',3);

set(tmph,'MarkerSize',30);
tmph=plot(interpvec(2,:),R0(lambsolve(:,2)'.^-1),'k-');    %bottom CI
set(tmph,'linewidth',3,'MarkerSize',30);
 


    %intersection with mock data
    tmph=plot(1./[real_lambda real_lambda], [R0(pt1^(-1)) R0(pt2^(-1))], 'r-');
    set(tmph,'linewidth',3);
    tmph=plot([0 1./real_lambda], [R0(pt2^(-1)) R0(pt2^(-1))], 'r--');
    set(tmph,'linewidth',3);
    tmph=plot([0 1./real_lambda], [R0(pt1^(-1)) R0(pt1^(-1))], 'r--');
    set(tmph,'linewidth',3);


hold off


%Plot Style
%title(strcat(titlename),'fontsize',20,'interpreter','latex');
%title('$\lambda^{-1}$ uncertainty mapped to $\mathcal{R}_0$','fontsize',20,'interpreter','latex');
xlabel('$\tau_c$ measured','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('$\mathcal{R}_{0}$ range','fontsize',20,'verticalalignment','bottom','interpreter','latex');
set(tmph,'linewidth',3);
set(gca,'fontsize',20);

%xlim([.9*min(interpvec(2,:)) 1.1*max(interpvec(1,:))])
xlim([.9*min(interpvec(2,:)) 33])
ylim([.92*min(R0(lambsolve(:,3)'.^-1)) 1.05*max(R0(lambsolve(:,3)'.^-1))])

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