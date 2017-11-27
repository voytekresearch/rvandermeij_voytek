function rmr_abcaudecog_getpac



% get details
info = rmr_abcaudecog_info;
%info.datapath = '/Users/roemer/Work/Data/tmpdata/BidetCaulet_ECoG/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/abcaudecog/';


% sub select for amua differences between conditions of some sort
%info.subj = {'GP15','GP22','GP28','GP35'};

% set
centerampfreq = 'yes';

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  % check for already done
  fn = [];
  fn{1} = [info.savepath currsubj '_' 'wplf_9to60'           '_' 'centerae' centerampfreq  '.mat'];
  fn{2} = [info.savepath currsubj '_' 'wplf_9to60to70to130'  '_' 'centerae' centerampfreq  '.mat'];
  fn{3} = [info.savepath currsubj '_' 'wplf_9to60toHFA'      '_' 'centerae' centerampfreq  '.mat'];
  %   if all(cellfun(@exist,fn,repmat({'file'},[1 numel(fn)])))
  %     continue
  %   end
  
  % fetch and preprocess data
  cfg = [];
  cfg.demean      = 'yes';
  cfg.prestim     = 0;
  cfg.poststim    = 0;
  cfg.reref       = 'yes';
  cfg.refchannel  = 'all';
  data = rmr_abcaudecog_fetchpreprocessdata(cfg,currsubj,info);
  
  
  % set trial bool
  % lo pitch
  attRRlostd = data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1; % attend right, right stim, standard,  correctly rejected, lo pitch
  attRLlostd = data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1; % attend right, left  stim, standard,  correctly rejected, lo pitch
  attLRlostd = data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1; % attend left,  right stim, standard,  correctly rejected, lo pitch
  attLLlostd = data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1; % attend left,  left stim,  standard,  correctly rejected, lo pitch
  attBRlostd = data.trialinfo(:,1) == 3 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1; % attend bin,   right stim, standard,  correctly rejected, lo pitch
  attBLlostd = data.trialinfo(:,1) == 3 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1; % attend bin,   left stim,  standard,  correctly rejected, lo pitch
  % hi pitch
  attRRhistd = data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2; % attend right, right stim, standard,  correctly rejected, hi pitch
  attRLhistd = data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2; % attend right, left  stim, standard,  correctly rejected, hi pitch
  attLRhistd = data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2; % attend left,  right stim, standard,  correctly rejected, hi pitch
  attLLhistd = data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2; % attend left,  left stim,  standard,  correctly rejected, hi pitch
  attBRhistd = data.trialinfo(:,1) == 3 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2; % attend bin,   right stim, standard,  correctly rejected, hi pitch
  attBLhistd = data.trialinfo(:,1) == 3 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2; % attend bin,   left stim,  standard,  correctly rejected, hi pitch
  %
  allstdcr   = data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4; % standard,  correctly rejected
  
  
  
  % cut out data
  cfg = [];
  cfg.toilim = [0 0.4];
  datasel = ft_redefinetrial(cfg,data);
  cfg = [];
  cfg.demean  = 'yes';
  cfg.detrend = 'yes';
  cfg.trials  = find(allstdcr & ((cellfun(@numel,data.time).'./datasel.fsample)>0.380)); % lame
  datasel = ft_preprocessing(cfg,datasel);
  
  
  % prep PAC settings
  ampfreq        = {[9:1:20 24:4:60]  ,70:4:130           ,100:10:200                    };
  phasfreq       = {ampfreq{1}        ,ampfreq{1}         ,ampfreq{1}                    };
  ampfreqtimwin  = {3./ampfreq{1}     ,3./ampfreq{2}      ,ones(size(ampfreq{3})) .* .050};
  phasfreqtimwin = {3./phasfreq{1}    ,3./phasfreq{2}     ,3./phasfreq{3}};
  amptstohfats   = {'no'              ,'no'               ,'yes'};
  
  % all PACs
  for iset = 1:3
    if ~exist(fn{iset},'file')
      
      % get pac
      cfg = [];
      cfg.ampfreq         = ampfreq{iset};
      cfg.ampfreqtimwin   = ampfreqtimwin{iset};
      cfg.phasfreq        = phasfreq{iset};
      cfg.phasfreqtimwin  = phasfreqtimwin{iset};
      cfg.amptstohfats    = amptstohfats{iset};
      cfg.centerampfreq   = centerampfreq;
      cfg.phasfreqnoamp   = 'no';
      cfg.handlespecedges = 'liberal';
      cfg.normalization   = 'withintrials';
      cfg.keeptrials      = 'no';
      wplfdata = rmr_phaseamplitudecoupling(cfg,datasel);
      
      % save
      save(fn{iset},'wplfdata','-v7.3');
      clear wplfdata
    end
    
    % do the two partitions
    for ish = 1:2
      currfn = [fn{iset}(1:end-4) '_' 'shpart'  num2str(ish) '.mat'];
      if ~exist(currfn,'file')
        
        % get pac
        cfg = [];
        cfg.trials          = ish:2:numel(datasel.trial);
        cfg.ampfreq         = ampfreq{iset};
        cfg.ampfreqtimwin   = ampfreqtimwin{iset};
        cfg.phasfreq        = phasfreq{iset};
        cfg.phasfreqtimwin  = phasfreqtimwin{iset};
        cfg.amptstohfats    = amptstohfats{iset};
        cfg.centerampfreq   = centerampfreq;
        cfg.phasfreqnoamp   = 'no';
        cfg.handlespecedges = 'liberal';
        cfg.normalization   = 'withintrials';
        cfg.keeptrials      = 'no';
        wplfdata = rmr_phaseamplitudecoupling(cfg,datasel);
        
        % save
        save(currfn,'wplfdata','-v7.3');
        clear wplfdata
      end
    end % ish
    
  end % iset
end





















function playground




%%%%%%%%%%%%%%%%%
%%% BETWEEN
% fetch info
info = rmr_abcaudecog_info;


% design
centerampfreq = 'yes';

for     isubj = 1  :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  fn = [];
  fn{1} = [info.savepath currsubj '_' 'wplf_9to60'           '_' 'centerae' centerampfreq  '.mat'];
  fn{2} = [info.savepath currsubj '_' 'wplf_9to60to70to130'  '_' 'centerae' centerampfreq  '.mat'];
  fn{3} = [info.savepath currsubj '_' 'wplf_9to60toHFA'      '_' 'centerae' centerampfreq  '.mat'];

  % skip if needed
  if ~exist(fn{3},'file')
    continue
  end
  
  % for each of the freqs, plot the tfr
  for ipac = 1:numel(fn)
    
    % load that shit up
    load(fn{ipac})
    
    % get layout
    load([info.savepath currsubj '_fieldtrip_layout.mat'])
    
    % plot
    if numel(wplfdata.ampfreq)>2
      figname = [currsubj '-' 'oscosc'];
    else
      figname = [currsubj '-' 'oschfa'];
    end
    
    % browse
    cfg = [];
    cfg.layout       = lay;
    cfg.windowname   = figname;
    cfg.plotlabel    = 'yes';
    cfg.plotoutline  = 'no';
    rmr_browsewplfdata(cfg,wplfdata);
    
  end % iperiod
  
  
end % isubj
%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%
%%% WITHIN
% fetch info
info = rmr_abcaudecog_info;

%info.datapath = '/Users/roemer/Work/Data/tmpdata/BidetCaulet_ECoG/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/abcaudecog/';

% design
centerampfreq = 'yes';

for     isubj = 1  :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  fn = [];
  fn{1} = [info.savepath currsubj '_' 'wplf_9to60'           '_' 'centerae' centerampfreq  '.mat'];
  fn{2} = [info.savepath currsubj '_' 'wplf_9to60to70to130'  '_' 'centerae' centerampfreq  '.mat'];
  fn{3} = [info.savepath currsubj '_' 'wplf_9to60toHFA'      '_' 'centerae' centerampfreq  '.mat'];
  
  
  % skip if needed
  if ~exist(fn{3},'file')
    continue
  end
  
  % for each of the freqs, plot the tfr
  for ipac = 1:numel(fn)
    
    % load that shit up
    load(fn{ipac})
    
    
    % get layout
    load([info.savepath currsubj '_fieldtrip_layout.mat'])
    
    
    % plot
    if numel(wplfdata.ampfreq)>2
      figname = [currsubj '-' 'oscosc'];
    else
      figname = [currsubj '-' 'oschfa'];
    end
    figure('numbertitle','off','name',figname)
    hold on
    
    % get clim
    nchan = numel(wplfdata.label);
    clim = [0 0];
    for ichan = 1:nchan
      clim(2) = max(clim(2),max(max(abs(squeeze(wplfdata.wplf(ichan,ichan,:,:))))));
    end
    clim = round(clim*100)./100;
    %clim = [0 0.1];
    
    % flg
    hashfa = size(wplfdata.wplf,3)==1;
    
    % rearrange wplf
    wplfold   = wplfdata.wplf;
    nfreqamp  = size(wplfold,3);
    nfreqphas = size(wplfold,4);
    wplf = zeros(nchan,nfreqamp,nfreqphas);
    for ichan = 1:nchan
      wplf(ichan,:,:) = wplfold(ichan,ichan,:,:);
    end
    
    % fake freqdata
    fakedata = [];
    fakedata.powspctrm = squeeze(abs(wplf));
    if ~hashfa
      fakedata.time   = wplfdata.phasfreq;
      fakedata.freq   = wplfdata.ampfreq;
      fakedata.dimord = 'chan_freq_time';
    else
      fakedata.freq = wplfdata.phasfreq;
      fakedata.dimord = 'chan_freq';
    end
    fakedata.label = wplfdata.label;
    
    % plot
    if hashfa
      cfg = [];
      cfg.layout = lay;
      cfg.vlim   = clim;
      ft_multiplotER(cfg,fakedata);
    else
      cfg = [];
      cfg.layout = lay;
      cfg.zlim   = clim;
      ft_multiplotTFR(cfg,fakedata);
    end
    
    
    
    %       % prep ticks
    %       ampfreq   = wplfdata{iset}.ampfreq';
    %       phasfreq  = wplfdata{iset}.phasfreq';
    %       nampfreq  = numel(ampfreq);
    %       nphasfreq = numel(phasfreq);
    %       amptickind    = 1:floor(nampfreq/6):nampfreq;
    %       ampticklabel  = num2str(round(ampfreq(amptickind)));
    %       phastickind   = 1:floor(nphasfreq/6):nphasfreq;
    %       phasticklabel = num2str(round(phasfreq(phastickind)));
    %
    %       % draw data
    %       label = wplfdata{iset}.label;
    %       nchan = numel(label);
    %       for ichan = 1:nchan
    %         subplot_tight(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan,0.025)
    %         %
    %         C = abs(squeeze(wplfdata{iset}.wplf(ichan,ichan,:,:)));
    %         imagesc(C);
    %         axis xy
    %         caxis(clim)
    %         title(wplfdata{iset}.label{ichan})
    %         if ichan == nchan
    %           set(gca,'ytick',amptickind,'yticklabel',ampticklabel,'xtick',phastickind,'xticklabel',phasticklabel)
    %           %colorbar('location','eastoutside')
    %         else
    %           set(gca,'ytick',amptickind,'xtick',amptickind,'yticklabel',[],'xticklabel',[])
    %         end
    %       end
    %       subplot_tight(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan+1,0.025)
    %       caxis(clim)
    %       colorbar
    %       axis off
    
    
  end % ipac
  
  
end % isubj
%%%%%%%%%%%%%%



































