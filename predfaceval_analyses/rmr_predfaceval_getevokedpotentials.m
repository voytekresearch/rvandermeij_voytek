function rmr_predfaceval_getevokedpotentials





% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

for     isubj = 1:1:numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  disp(['working on ' currsubj])
  
  % for artifacts been detected
  if ~info.(currsubj).trlartfctflg
    disp(['artifacts not yet detected for ' currsubj])
    continue
  end
  % check for alreadt done
  fn = [info.savepath currsubj '_' 'timelockspoststim_001hplp10' '.mat'];
  if exist(fn,'file');
    continue
  end
  
  
  % combine data over sessions
  %   fndat = [info.savepath currsubj '_' 'allsess_bs540lp_nodetrend_500mspre_700mspost' '.mat'];
  %   if ~exist(fndat,'file')
  % fetch and preprocess data
  cfg = [];
  cfg.demean      = 'yes';
  cfg.detrend     = 'no';
  cfg.prestim     = 0.5;
  cfg.poststim    = 0.5+0.2;
  cfg.reref       = 'yes';
  cfg.refchannel  = 'all';
  cfg.bsfilter    = 'yes';
  % override low pass and set highpass
  hdr = ft_read_header([info.datapath currsubj '/' info.sessiondata.(currsubj){1}]);
  % set the low pass
  cfg.lpfilter      = 'yes';
  cfg.lpfreq        = 10;
  cfg.lpfilttype    = 'firws';
  cfg.lpfiltdir     = 'onepass-zerophase';
  cfg.lpfiltwintype = 'blackman';
  cfg.lpfiltord     = round(1*hdr.Fs);
  cfg.usefftfilt    = 'yes';
  % set the high pass
  cfg.hpfilter      = 'yes';
  cfg.hpfreq        = 0.01;
  cfg.hpfilttype    = 'firws';
  cfg.hpfiltdir     = 'onepass-zerophase';
  cfg.hpfiltwintype = 'blackman';
  cfg.hpfiltord     = round(11*hdr.Fs);
  data = rmr_predfaceval_fetchpreprocessdata(cfg,currsubj,info);
  
  %     % save dat
  %     save(fndat,'data','-v7.3');
  %   else
  %     load(fndat)
  %   end
  
  % get poststim and baseline 200ms before that
  cfg = [];
  cfg.toilim = [0.8 1.9];
  data = ft_redefinetrial(cfg,data);
  cfg = [];
  cfg.offset = -1.*data.fsample;
  data = ft_redefinetrial(cfg,data);
  cfg = [];
  cfg.demean         = 'yes';
  cfg.baselinewindow = [-0.2 0];
  cfg.trials         = data.trialinfo(:,5)<1.5 & data.trialinfo(:,3)==1; % select only trials of same duration and that had a correct response;
  data = ft_preprocessing(cfg,data);
  
  % get trialdef
  valencetc = data.trialinfo(:,2); % fearful (1) or neutral (2) face trial
  predtc    = data.trialinfo(:,1); % pred (1) or unpred (2) face trial
  
  % get timelocks
  % 1 = fear
  % 2 = neut
  % 3 = fearpred
  % 4 = fearunpred
  % 5 = neutpred
  % 6 = neutunpred
  timelock = [];
  cfg = [];
  cfg.keeptrials = 'yes';
  cfg.vartrllength = 2;
  cfg.trials       = valencetc==1;
  timelock{1} = ft_timelockanalysis(cfg,data);
  cfg.trials       = valencetc==2;
  timelock{2} = ft_timelockanalysis(cfg,data);
  cfg.trials       = valencetc==1 & predtc==1;
  timelock{3} = ft_timelockanalysis(cfg,data);
  cfg.trials       = valencetc==1 & predtc==2;
  timelock{4} = ft_timelockanalysis(cfg,data);
  cfg.trials       = valencetc==2 & predtc==1;
  timelock{5} = ft_timelockanalysis(cfg,data);
  cfg.trials       = valencetc==2 & predtc==2;
  timelock{6} = ft_timelockanalysis(cfg,data);
  % save
  save(fn,'timelock','-v7.3');
  
end













function playground















% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

for     isubj = 2  :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  disp(['working on ' currsubj])
  fn = [info.savepath currsubj '_' 'timelockspoststim' '.mat'];
  %fn = [info.savepath currsubj '_' 'timelockspoststim_001hplp10' '.mat'];
  
  % skip if needed
  if ~exist(fn,'file')
    continue
  end
  
  % load that shit up
  load(fn)
  
  % select channels
  ind = [];
  type = {'RAM','LAM','AM'};
  for itype = 1:numel(type)
    ind = [ind; find(strncmp(timelock{1}.label,type{itype},numel(type{itype})))];
  end
  for itl = 1:6
    timelock{itl}.label = timelock{itl}.label(ind);
    timelock{itl}.avg   = timelock{itl}.avg(ind,:);
    timelock{itl}.var   = timelock{itl}.var(ind,:);
    timelock{itl}.dof   = timelock{itl}.dof(ind,:);
  end
  
  % set some things
  nchan = numel(timelock{1}.label);
  
  % 1 = fear
  % 2 = neut
  % 3 = fearpred
  % 4 = fearunpred
  % 5 = neutpred
  % 6 = neutunpre
  
  % plot fearneut
  figure('numbertitle','off','name',[currsubj 'fearneut'])
  for ichan = 1:nchan
    subplot(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan)
    hold on
    %
    sel1 = 1;
    sel2 = 2;
    %
    avg1 = timelock{sel1}.avg(ichan,:);
    avg2 = timelock{sel2}.avg(ichan,:);
    sem1 = sqrt(timelock{sel1}.var(ichan,:)) ./ sqrt(mean(timelock{sel1}.dof(:))); % dit gaat ervan uit dat alle trials, pre-averaging, even lang waren
    sem2 = sqrt(timelock{sel2}.var(ichan,:)) ./ sqrt(mean(timelock{sel2}.dof(:)));
    time = timelock{sel1}.time;
    % plot patch, then lines
    x = [time time(end:-1:1)];
    y = [avg1-sem1 avg1(end:-1:1)+sem1(end:-1:1)];
    patch(x,y,'b','edgecolor','none','facealpha',0.15)  % patch wil een set x en y coordinaten hebben, die de 'randen' van het object aangeven
    y = [avg2-sem2 avg2(end:-1:1)+sem2(end:-1:1)];
    patch(x,y,'gr','edgecolor','none','facealpha',0.15)  % een hoger getal hier (max = 1) maakt de sems minder transparant
    plot(time,avg1,'b') % door het gemiddelde zelf pas later te plotten, staat hij niet achter de patch, maar op de voorgrond
    plot(time,avg2,'gr')
    title(timelock{1}.label{ichan})
    % plot axis
    xlim = [-0.2 0.9]; % x-as limit
    ylim = [-60 60]; %
    set(gca,'xlim',xlim,'ylim',ylim)
    line([0 0],ylim,'color','k') % deze twee regels plotten de x en y as door t=0 en amp=0 heen
    line(xlim,[0 0],'color','k')
    xlabel('time');ylabel('voltage')
    set(gca,'ytick',[ylim(1) 0 ylim(2)],'xtick',[xlim(1) 0 xlim(2)]) % dit zorgt ervoor dat de assen 'streepjes', de ticks, beperkt zijn
  end
  
  
  
end % isubj
















timelock

cfg = [];
cfg.layout = 'ordered';
lay = ft_prepare_layout(cfg,timelock);

% 'fake' a freq structure with powspctrm from the timelock structure (created with cfg.keeptrials = 'yes')
freq = timelock;
freq.freq = 1:size(freq.trial,1);
freq.powspctrm = permute(freq.trial,[2 1 3]); % dimension order of 'trial' field is different than powspctrm
freq.dimord = 'chan_freq_time';
% (leave other fields in the data, they'll be left alone)

figure
cfg = [];
cfg.layout      = '';
cfg.parameter   = 'powspctrm';
cfg.interactive = 'yes';
cfg.zlim        = 'maxabs'; % colorscale: -max of the abs to 0 to max of the abs
ft_multiplotTFR(cfg,freq)



close all































