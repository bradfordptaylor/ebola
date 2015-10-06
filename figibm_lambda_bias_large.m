clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figibm_lambda_bias_large';
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

% main data goes here
info.n_E    = 1;        % Number of stages for Gamma
info.T_E    = 11;       % Days exposed, on average
info.b_E    = info.n_E/info.T_E;  % Control parameter for Gamma
info.T_I    = 6;        % Number of days infectious
info.b_I    = 1/info.T_I;% Control parameter for I
info.T_D    = 3;        % Number of days until burial
info.b_D    = 1./info.T_D;% Control parameter for D
info.fracR0_dead_range = 0:0.1:1;  % Range for the evaluation
info.sigma=1/info.T_E;
info.gamma=1/info.T_I;
info.f=0.7;
info.N=10^6;
info.Itrigger=50;
info.rho=info.b_D;
info.num_ensembles=1000;
info.tmax=500;
info.t0=0;

% Now perform the search
info.lambda_target=1./[14 21 28];

fp = fopen('figibm_lambda_bias_large.mat');
if (fp<0)
  %
  % Outer loop - lambda
  %
  more off
  for li = 1:length(info.lambda_target),
    li
    info.lambda = info.lambda_target(li);
    z=-info.lambda;
    % Shared by all in the next loop
    M_E = (info.b_E/(info.b_E-z))^info.n_E;
    M_I = info.b_I/(info.b_I-z);
    M_D = info.b_D/(info.b_D-z);
    info.M_E=M_E;
    info.M_I=M_I;
    info.M_D=M_D;
  
    %
    % Intermediate loop - fraction of transmission due to burials
    %
    for rhoi = 1:length(info.fracR0_dead_range),
      rhoi
      info.fracR0_dead = info.fracR0_dead_range(rhoi);
      info.c_I = 1-info.fracR0_dead;
      info.Mtot = info.c_I*M_E*M_I+info.fracR0_dead*M_E*M_I*M_D;
      info.R0 = 1/info.Mtot;
      stats(li,rhoi).R0=info.R0;
      stats(li,rhoi).fracR0_dead = info.fracR0_dead;
      stats(li,rhoi).betaD = info.fracR0_dead.*info.R0/info.f/info.T_D;
      stats(li,rhoi).betaI = (1-info.fracR0_dead).*info.R0/info.T_I;
      stats(li,rhoi).R0_seir= 1/(M_E*M_I);
      stats(li,rhoi).beta_seir = stats(li,rhoi).R0_seir/info.T_I;
      info.beta_D=stats(li,rhoi).betaD;
      info.beta_I=stats(li,rhoi).betaI;
      stats(li,rhoi).info=info;
    
      %
      % Final loop - ensembles of trajectories
      %
      for n = 1:info.num_ensembles,
        n
        clear gstats
        t = 0;
        % S E I R D B
        y = [info.N-1 0 1 0 0 0]; 
        gstats.t=t;
        gstats.y(1,:)=y;
        nume = 1;
        cumI = 0;
        tmax=info.tmax;
        ttrigger=0;
        while (t<=tmax & sum(y(1:3))>0 & sum(y([2 3 5])>0))  % Epidemic is still possible
          e_list = calc_elist(y,info);
          t_next = t-1/sum(e_list)*log(rand);
          y_new = perform_event(y,info,e_list);
          nume = nume+1;
          gstats.t(nume)=t_next;
          gstats.y(nume,:)=y_new;
          cumI = y(3);
          if (cumI>info.Itrigger & ttrigger==0)
            tmax=t+2.5/info.lambda;
            ttrigger=1;
          end
          y = y_new;
          t = t_next;
        end
        
        % Do the sampling, daily
        gstats.cens_dt = 1;
        gstats.cens_t = info.t0:gstats.cens_dt:info.tmax;
        tmpfloort = floor(gstats.t);
        for i=1:length(gstats.cens_t),
          tmpi = find(tmpfloort == gstats.cens_t(i));
          if (~isempty(tmpi))
            if length(tmpi)>1
              gstats.cens_y(i,:) = floor(median(gstats.y(tmpi,:)));
            else
              gstats.cens_y(i,:) = gstats.y(tmpi,:);
            end
          else
            tmpi = find(gstats.t < gstats.cens_t(i));
            gstats.cens_y(i,:) = gstats.y(tmpi(end),:);
          end
        end
  
        % Do the fitting
        tmpi = find (gstats.cens_y(:,3)>info.Itrigger);
        if (~isempty(tmpi))
          tmpj = find (gstats.cens_t<=tmpi(1)+2/info.lambda & gstats.cens_t>=gstats.cens_t(tmpi(1)));
          [p,s]=polyfit(gstats.cens_t(tmpj)',log(gstats.cens_y(tmpj,3)),1);
          stats(li,rhoi).lambda(n) =p(1);
          stats(li,rhoi).offset(n) = p(2);
        else
          stats(li,rhoi).lambda(n) =-1;
          stats(li,rhoi).offset(n) =-1;
        end
  
      end  % Final loop
  
    end    % Intermediate loop
   
  end	 % Outer loop
  
  % Save before plotting
  save figibm_lambda_bias_large stats

else
  load('figibm_lambda_bias_large.mat');
end

% Now plot
info=stats(1,1).info;
for li = 1:length(info.lambda_target),
  % Next loop
  for rhoi = 1:length(info.fracR0_dead_range),
    info = stats(li,rhoi).info;
    lambda = info.lambda;
    tmpi=find(stats(li,rhoi).lambda>0);   
    numl = length(tmpi);
    llow = floor(numl*0.025);
    lhigh = floor(numl*0.975);
    lsort = sort(stats(li,rhoi).lambda(tmpi),'ascend');
    stats(li,rhoi).lambda_low = lsort(llow); 
    stats(li,rhoi).lambda_high = lsort(lhigh); 
    mlambda(rhoi) = mean(stats(li,rhoi).lambda(tmpi));
    stdlambda(rhoi) = std(stats(li,rhoi).lambda(tmpi));
    mchartime(rhoi) = mean(1./stats(li,rhoi).lambda(tmpi));
    stdchartime(rhoi) = std(1./stats(li,rhoi).lambda(tmpi));
    stats(li,rhoi).lambda_mean=mlambda(rhoi);
    stats(li,rhoi).mchartime=mchartime(rhoi);
    % Plot the middle
    tmph=plot(info.fracR0_dead_range(rhoi),mchartime(rhoi),'ko');
    set(tmph,'markerfacecolor',[0.0 0.0 0.0]+(li-1)*[0.3 0.3 0.3],'markersize',12);
    hold on
    % Then plot the caps
    tmph=plot(info.fracR0_dead_range(rhoi),1/stats(li,rhoi).lambda_low,'kv');
    set(tmph,'markerfacecolor',[0.0 0.0 0.0]+(li-1)*[0.3 0.3 0.3],'markersize',12);
    % Then plot the caps
    tmph=plot(info.fracR0_dead_range(rhoi),1/stats(li,rhoi).lambda_high,'k^');
    set(tmph,'markerfacecolor',[0.0 0.0 0.0]+(li-1)*[0.3 0.3 0.3],'markersize',12);
    % Now connect
    tmph=plot(info.fracR0_dead_range(rhoi)*ones(1,3),[1/stats(li,rhoi).lambda_low mchartime(rhoi) 1/stats(li,rhoi).lambda_high],'k--');
    set(tmph,'linewidth',2,'color',[0.0 0.0 0.0]+(li-1)*[0.3 0.3 0.3]);
  end
  % Draw the ideal
  tmph=plot(info.fracR0_dead_range,(1/stats(li,1).info.lambda)*ones(1,length(info.fracR0_dead_range)),'k-');
  hold on
  set(tmph,'linewidth',3);
  set(tmph,'color',[0.0 0.0 0.0]+(li-1)*[0.2 0.2 0.2]);
end

% loglog(,, '');
%
%
% Some helpful plot commands
% tmph=plot(x,y,'ko');
% set(tmph,'markersize',10,'markerfacecolor,'k');
% tmph=plot(x,y,'k-');
% set(tmph,'linewidth',2);

set(gca,'fontsize',20);

% for use with layered plots
% set(gca,'box','off')

% adjust limits
% tmpv = axis;
% axis([]);
% ylim([]);
xlim([-0.05 1.05]);
ylim([10 40]);

% change axis line width (default is 0.5)
% set(tmpa1,'linewidth',2)

% fix up tickmarks
% set(gca,'xtick',[1 100 10^4])
% set(gca,'ytick',[1 100 10^4])

% creation of postscript for papers
% psprint(tmpxfigfilename);

% the following will usually not be printed 
% in good copy for papers
% (except for legend without labels)

% legend
% tmplh = legend('stuff',...);
% tmplh = legend('','','');
% remove box
% set(tmplh,'visible','off')
% legend('boxoff');

xlabel('$\rho_D={\cal{R}}_0(\mathrm{dead~})/{\cal{R}}_0$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Characteristic time, $\langle 1/\lambda\rangle$ (days)','fontsize',20,'verticalalignment','bottom','interpreter','latex');
% title('','fontsize',24)
% 'horizontalalignment','left');

% for writing over the top
% coordinates are normalized again to (0,1.0)
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
