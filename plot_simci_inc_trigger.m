%Plot country results -- Guinea
close all
clear all

%load processed data
filenameload = strcat('data/ci_sim_inc_trig50.mat');
load(filenameload)

%solve for "measured" rate constant

real_lambda = 1/20;  % mock data


%Find full intersection of measured rate constant with CI from sims.
%Intersection estimated by linear interpolation


    aboveptlow = find(diff(sign(1./data.poislambda_low-1./real_lambda)));
    abovepthigh = find(diff(sign(1./data.poislambda_high-1./real_lambda)));
    yone = [data.poislambda_low(aboveptlow) data.poislambda_high(abovepthigh)].^(-1);
    ytwo = [data.poislambda_low(aboveptlow+1) data.poislambda_high(abovepthigh+1)].^(-1);
    xone = [data.xsupport(aboveptlow) data.xsupport(abovepthigh)].^(-1);
    xtwo = [data.xsupport(aboveptlow+1) data.xsupport(abovepthigh+1)].^(-1);
    xsolve = ((1./real_lambda-yone).*(xtwo-xone)./(ytwo-yone))+xone;
    
    savename = strcat('data/xsolve_simci_inc');
    save(savename,'xsolve')



%lattice method in obtaining simulated CI for each measured rate constant
latticenum = 500;
%set and save x-axis for R0 plot, x-axis for upper CI and lower CI differ
interpvec(1,:) = linspace(min(1./data.poislambda_low),max(1./data.poislambda_low),latticenum); %measured for high CI in lambda-lambda plot, this becomes support
interpvec(2,:) = linspace(min(1./data.poislambda_high),max(1./data.poislambda_high),latticenum); %measured for low CI in lambda-lambda plot, this becomes support
interpvec(3,:) = linspace(min(1./data.poismlambda),max(1./data.poismlambda),latticenum);
savename = strcat('data/interpvec_simci_inc');
save(savename,'interpvec')

%loop over lattice number (different measured rate constant for upper and
%lower CI) and find boundary value

for jj=1:length(interpvec)
    pt1 = interpvec(1,jj); %High lambda CI support point
    pt2 = interpvec(2,jj); %low lambda CI support point
    pt3 = interpvec(3,jj);
    
    %See comments in figure 3 for explanation of code
    [checkbot,locbot] = ismember(pt1,1./data.poislambda_low);
    if checkbot>0
        aboveptlow = locbot;
    else
        aboveptlow = find(diff(sign(1./data.poislambda_low-pt1)));
        aboveptlow = max(aboveptlow);
    end
     
    [checktop,loctop] = ismember(pt2,1./data.poislambda_high);
    if sum(checktop)>0
        abovepthigh = loctop;
    else
        abovepthigh = find(diff(sign(1./data.poislambda_high-pt2)));
        abovepthigh = max(abovepthigh);
    end
    
    [checkmed,locmed] = ismember(pt3,1./data.poismlambda);
    if sum(checkmed)>0
        aboveptmed = locmed;
    else
        aboveptmed = find(diff(sign(1./data.poismlambda-pt3)));
        aboveptmed = max(aboveptmed);
    end

    if aboveptlow~=length(data.poislambda_low)
        yone = data.poislambda_low(aboveptlow).^(-1);
        ytwo = data.poislambda_low(aboveptlow+1).^(-1);
        xone = data.xsupport(aboveptlow).^(-1);
        xtwo = data.xsupport(aboveptlow+1).^(-1);
    else
        yone = data.poislambda_low(aboveptlow-1).^(-1);
        ytwo = data.poislambda_low(aboveptlow).^(-1);
        xone = data.xsupport(aboveptlow-1).^(-1);
        xtwo = data.xsupport(aboveptlow).^(-1);
    end
    
    if abovepthigh~=length(data.poislambda_high)
        yone(2) = data.poislambda_high(abovepthigh).^(-1);
        ytwo(2) = data.poislambda_high(abovepthigh+1).^(-1);
        xone(2) = data.xsupport(abovepthigh).^(-1);
        xtwo(2) = data.xsupport(abovepthigh+1).^(-1);
    else
        yone(2) = data.poislambda_high(abovepthigh-1).^(-1);
        ytwo(2) = data.poislambda_high(abovepthigh).^(-1);
        xone(2) = (data.xsupport(abovepthigh-1)).^(-1);
        xtwo(2) = (data.xsupport(abovepthigh)).^(-1);
         
    end
    
    if aboveptmed~=length(data.poismlambda)
        yone(3) = data.poismlambda(aboveptmed).^(-1);
        ytwo(3) = data.poismlambda(aboveptmed+1).^(-1);
        xone(3) = data.xsupport(aboveptmed).^(-1);
        xtwo(3) = data.xsupport(aboveptmed+1).^(-1);
    else
        yone(3) = data.poismlambda(aboveptmed-1).^(-1);
        ytwo(3) = data.poismlambda(aboveptmed).^(-1);
        xone(3) = (data.xsupport(aboveptmed-1)).^(-1);
        xtwo(3) = (data.xsupport(aboveptmed)).^(-1);
         
    end
    lambsolve(jj,:) = (([pt1 pt2 pt3]-yone).*(xtwo-xone)./(ytwo-yone))+xone;
end

savename = strcat('data/lambsolve_simci_inc');
save(savename,'lambsolve')
 

clf;
% automatically create postscript whenever
% figure is drawn
eval(strcat('tmpfilename = ''figures/fig_simci_inc_large'';'))
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
patch([1./data.xsupport fliplr(1./data.xsupport)],[1./data.poislambda_low fliplr(1./data.poislambda_high)],[.8 .8 .8])

hold on
tmph = plot(1./data.xsupport,1./data.poismlambda,'b.-');
set(tmph,'linewidth',3,'MarkerSize',30);

 tmph=plot(1./data.xsupport,1./data.poislambda_low,'k.-');
 set(tmph,'linewidth',3,'MarkerSize',30);
 tmph=plot(1./data.xsupport,1./data.poislambda_high,'k.-');
 set(tmph,'linewidth',3,'MarkerSize',30);

 tmph=plot(xsolve, [1./real_lambda 1./real_lambda], 'r-');
  set(tmph,'linewidth',3);
  tmph=plot([xsolve(1) xsolve(1)], [0 1./real_lambda], 'r--');
  set(tmph,'linewidth',3);
  tmph=plot([xsolve(2) xsolve(2)], [0 1./real_lambda], 'r--');
  set(tmph,'linewidth',3);
  


%Plot style
ylim([min(1./data.poislambda_high) 33]);%1.01*max(1./data.lambda_low)])
xlim(ylim);	%axis limits to show 1-1 line
%  set(tmph,'linewidth',3);
%    tmph=plot(xlim, [1./real_lambda 1./real_lambda], 'r-');
%    set(tmph,'linewidth',3);
    hold off
%title(strcat(titlename),'fontsize',20,'interpreter','latex');
xlabel('$\tau_c$ theory','fontsize',20,'verticalalignment','top','interpreter','latex');
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