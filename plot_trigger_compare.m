function plot_trigger_compare


datastring = 'cum';
allxsolve = zeros(2,5);
for mm=2:5
    Imin = mm*10;
    fileloadcum = strcat('data/ci_sim_',datastring,'_trig',num2str(Imin));
    load(fileloadcum)





%solve for "measured" rate constant

real_lambda = 1/20;  %mock data






%Find full intersection of measured rate constant with CI for sierra leone
%and liberia. For guinea, find intersection with upper CI. Same methodology
%as in figure 3 and figure 4.


    aboveptlow = find(diff(sign(1./data.poislambda_low-1./real_lambda)));
    abovepthigh = find(diff(sign(1./data.poislambda_high-1./real_lambda)));
    yone = [data.poislambda_low(aboveptlow) data.poislambda_high(abovepthigh)].^(-1);
    ytwo = [data.poislambda_low(aboveptlow+1) data.poislambda_high(abovepthigh+1)].^(-1);
    xone = [data.xsupport(aboveptlow) data.xsupport(abovepthigh)].^(-1);
    xtwo = [data.xsupport(aboveptlow+1) data.xsupport(abovepthigh+1)].^(-1);
    xsolve = ((1./real_lambda-yone).*(xtwo-xone)./(ytwo-yone))+xone;
    allxsolve(:,mm) = xsolve;
    %savename = strcat('data/xsolve_simci_cum');
    %save(savename,'xsolve')



end

%savename = strcat('data/lambsolve_simci_cum');
%save(savename,'lambsolve')
 

clf;
% automatically create postscript whenever
% figure is drawn
eval(strcat('tmpfilename = ''figures/fig_triggercompare'';'))
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
for pp = 2:5
tmph = plot([pp pp]*10,allxsolve(:,pp),'b-');
set(tmph,'linewidth',3,'MarkerSize',30);
hold on
tmph = plot([1 -1]+[pp pp]*10,allxsolve(1,pp).*ones(1,2),'b-');
set(tmph,'linewidth',3);
tmph = plot([1 -1]+[pp pp]*10,allxsolve(2,pp)*ones(1,2),'b-');
set(tmph,'linewidth',3);
end
hold off
xlim([10 60])


%Plot style

xlabel('Cumulative trigger population','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('$\tau_c$ CI','fontsize',20,'verticalalignment','bottom','interpreter','latex');
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