function rmr_abcaudecog_gethfa



% get details
info = rmr_abcaudecog_info;

% set stuff
savepath   = [info.savepath 'hfa/'];

% sub select, get rid of the JH's
info.subj = {'GP15','GP22','GP28','GP35','ST1','ST6','ST8'};

% sub select for hfa differences between conditions of some sort
%info.subj = {'GP15','GP22','GP28','GP35'};

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  % check for already done
  if numel(info.(currsubj).lohipitch)==2
    fn1 = [savepath currsubj '_' 'HFA100-10-200_attign' '_' info.(currsubj).lohipitch{1} '-pitch' '.mat'];
    fn2 = [savepath currsubj '_' 'HFA100-10-200_attign' '_' info.(currsubj).lohipitch{2} '-pitch' '.mat'];
    fn3 = [savepath currsubj '_' 'HFA100-10-200_pitchbin' '.mat'];
    fn4 = [savepath currsubj '_' 'HFA100-10-200_pitchmaxtrials' '.mat'];
    if all([exist(fn1,'file') exist(fn2,'file') exist(fn3,'file') exist(fn4,'file')])
      continue
    end
  else
    fn1 = [savepath currsubj '_' 'HFA100-10-200_attign' '_' info.(currsubj).lohipitch{1} '-pitch' '.mat'];
    fn2 = [savepath currsubj '_' 'HFA100-10-200_pitchbin' '.mat'];
    fn3 = [savepath currsubj '_' 'HFA100-10-200_pitchmaxtrials' '.mat'];
    if all([exist(fn1,'file') exist(fn2,'file') exist(fn3,'file')])
      continue
    end
  end
  
  
  datacmb = [];
  for    isess = 1   :numel(info.session.(currsubj))
    % set
    currsess = info.session.(currsubj){isess};
    
    % fetch data
    cfg = [];
    cfg.datapath = [info.datapath currsubj '/data/' currsess '/'];
    cfg.subject  = currsubj;
    cfg.session  = currsess;
    data = rmr_abcaudecog_data2fieldtrip(cfg);
    % apply filters
    cfg = [];
    cfg.bsfilter   = 'no';
    cfg.bsfilttype = 'but';
    cfg.bsfiltord  = 2;
    bsfiltpeak  = info.(currsubj).bsfilt.peak;
    bsfiltbandw = info.(currsubj).bsfilt.halfbandw;
    cfg.bsfreq     = [bsfiltpeak-bsfiltbandw;  bsfiltpeak+bsfiltbandw]';
    if data.fsample==1000
      cfg.lpfilter   = 'yes';
      cfg.lpfreq     = 350;
      cfg.lpfilttype = 'fir';
      cfg.lpfiltord  = 1000;
    else
      cfg.lpfilter   = 'yes';
      cfg.lpfreq     = 540;
      cfg.lpfilttype = 'fir';
      cfg.lpfiltord  = 1000;
    end
    %       cfg.lpfreq     = 250;
    %       cfg.lpfilttype = 'but';
    %       cfg.lpfiltord  = 25; % stable at 3000 and 1000 samping fs
    cfg.demean     = 'yes'; % there's a huge offset, get rid of it
    cfg.reref      = 'yes';
    cfg.refchannel = 'all';
    cfg.channel    = setdiff(data.label,info.(currsubj).badchan);
    data = ft_preprocessing(cfg,data);
    % replace labels
    cfg = [];
    cfg.subject = currsubj;
    data = rmr_abcaudecog_replacechanlabels(cfg,data);
    % define trials
    cfg = [];
    cfg.subj    = currsubj;
    cfg.isess   = isess;
    cfg.fsample = data.fsample;
    trl = rmr_abcaudecog_definetrialsfrompos(cfg);
    % add 200ms baseline
    trl(:,1) = trl(:,1) - round(0.2 .* data.fsample);
    trl(:,3) = -round(0.2 .* data.fsample);
    % cutout trials because of filtering
    remind = find((trl(:,1)./data.fsample)<info.(currsubj).bsfilt.edgeartlen);
    remind = [remind; find((trl(:,2)./data.fsample)>(data.time{1}(end)-info.(currsubj).bsfilt.edgeartlen))];
    trl(remind,:) = [];
    % redefine trials in data structure
    cfg = [];
    cfg.trl = trl;
    data = ft_redefinetrial(cfg,data);
    % rejects artifact
    load([info.savepath 'artfctdef/' currsubj '_' currsess '_artfctdef.mat'])
    cfg = [];
    cfg.artfctdef = artfctdef;
    data = ft_rejectartifact(cfg,data);
    % keep data for combining
    data = rmfield(data,'sampleinfo');
    datacmb{isess} = data;
  end
  % combine sessions
  if numel(info.session.(currsubj))>1
    data = ft_appenddata([],datacmb{:});
    clear datacmb
  else
    data = datacmb{1};
    clear datacmb
  end
  
  
  % set trial indices
  % lo pitch
  attRRlostd = find(data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1); % attend right, right stim, standard,  correctly rejected, lo pitch
  attRLlostd = find(data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1); % attend right, left  stim, standard,  correctly rejected, lo pitch
  attLRlostd = find(data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1); % attend left,  right stim, standard,  correctly rejected, lo pitch
  attLLlostd = find(data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1); % attend left,  left stim,  standard,  correctly rejected, lo pitch
  attBRlostd = find(data.trialinfo(:,1) == 3 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1); % attend bin,   right stim, standard,  correctly rejected, lo pitch
  attBLlostd = find(data.trialinfo(:,1) == 3 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 1); % attend bin,   left stim,  standard,  correctly rejected, lo pitch
  % hi pitch
  attRRhistd = find(data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2); % attend right, right stim, standard,  correctly rejected, hi pitch
  attRLhistd = find(data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2); % attend right, left  stim, standard,  correctly rejected, hi pitch
  attLRhistd = find(data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2); % attend left,  right stim, standard,  correctly rejected, hi pitch
  attLLhistd = find(data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2); % attend left,  left stim,  standard,  correctly rejected, hi pitch
  attBRhistd = find(data.trialinfo(:,1) == 3 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2); % attend bin,   right stim, standard,  correctly rejected, hi pitch
  attBLhistd = find(data.trialinfo(:,1) == 3 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4 & data.trialinfo(:,7) == 2); % attend bin,   left stim,  standard,  correctly rejected, hi pitch

  
  % set toi sampling
  [dum maxind] = max(cellfun(@numel,data.time));
  toi = data.time{maxind}(1:1:end);
  %     if data.fsample>3000
  %       [dum maxind] = max(cellfun(@numel,data.time));
  %       toi = data.time{maxind}(1:8:end);
  %     else
  %       [dum maxind] = max(cellfun(@numel,data.time));
  %       toi = data.time{maxind}(1:2:end);
  %     end
  
  % get filtdata via freqanalysis with 100ms windows, per freq
  freqoi = 100:10:200;
  nchan  = numel(data.label);
  ntrial = numel(data.trial);
  ntime  = numel(toi);
  nfreq  = numel(freqoi);
  filtdata = data;
  for itrial = 1:ntrial
    filtdata.trial{itrial} = zeros(nchan,ntime);
    filtdata.time{itrial}  = toi;
  end
  for ifreq = 1:nfreq
    disp(['getting HFA ampspctrm of ' num2str(freqoi(ifreq)) 'Hz'])
    % get powspctrm
    cfg = [];
    cfg.feedback   = 'none';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.pad        = 1; % pad out to 1 second;
    cfg.method     = 'mtmconvol';
    cfg.output     = 'pow';
    cfg.toi        = toi;
    cfg.foi        = freqoi(ifreq);
    cfg.taper      = 'hanning';
    cfg.t_ftimwin  = ones(numel(cfg.foi),1) .* .100; % THE BELOW ASSUMES identical time-windows for each frequency
    freqdata = ft_freqanalysis(cfg,data);
    ampenvspctrm = squeeze(sqrt(freqdata.powspctrm));  
    % zscore and merge as a running sum
    tmp = reshape(permute(ampenvspctrm(:,:,:),[2 1 4 3]),[nchan ntrial*ntime]);
    meanampenv = repmat(squeeze(nanmean(tmp,2)),[1 ntime]);
    stdampenv  = repmat(squeeze(nanstd(tmp,0,2)),[1 ntime]);
    for itrial = 1:ntrial
      currfreqtrial = permute(ampenvspctrm(itrial,:,:),[2 3 1]);
      currfreqtrial = (currfreqtrial-meanampenv)./stdampenv;  
      filtdata.trial{itrial} = filtdata.trial{itrial} + (currfreqtrial ./ numel(freqoi));
    end
  end
  % remove NaNs from filtdata (raw data cannot contain NaNs)
  for itrial = 1:ntrial
    nanind = isnan(sum(filtdata.trial{itrial},1));
    filtdata.trial{itrial} = filtdata.trial{itrial}(:,~nanind);
    filtdata.time{itrial}  = filtdata.time{itrial}(:,~nanind);
  end
  clear freqdata ampenvspctrm data tmp
  
  % baseline
  cfg = [];
  cfg.demean = 'yes';
  cfg.baselinewindow = [-0.1 0];
  filtdata = ft_preprocessing(cfg,filtdata);
  
  % set global timelocking settings
  cfgtl = [];
  cfgtl.keeptrials       = 'yes';
  cfgtl.covariance       = 'yes';
  cfgtl.covariancewindow = [0 toi(end)];
  cfgtl.vartrllength     = 2;
  
  %%%%% get tl for attention
  for iset = 1:numel(info.(currsubj).lohipitch)
    % determine which stimuli were (un)attended
    switch info.(currsubj).leftrightelec
      case 'right'
        if      strcmp(info.(currsubj).lohipitch{iset},'low')
          atttrials = attLLlostd;
          igntrials = attRLlostd;
        elseif  strcmp(info.(currsubj).lohipitch{iset},'high')
          atttrials = attLLhistd;
          igntrials = attRLhistd;
        else
          error('woopsie')
        end
        
      case 'left'
        if     strcmp(info.(currsubj).lohipitch{iset},'low')
          atttrials = attRRlostd;
          igntrials = attLRlostd;
        elseif strcmp(info.(currsubj).lohipitch{iset},'high')
          atttrials = attRRhistd;
          igntrials = attLRhistd;
        else
          error('woopsie')
        end
        
      otherwise
        error('woopsie')
    end
       
    % get timelocked average
    cfgtl.trials           = atttrials;
    tlattend = ft_timelockanalysis(cfgtl,filtdata);
    cfgtl.trials           = igntrials;
    tlignore = ft_timelockanalysis(cfgtl,filtdata);
    
    % save shit
    fn = [savepath currsubj '_' 'HFA100-10-200_attign' '_' info.(currsubj).lohipitch{iset} '-pitch' '.mat'];
    save(fn,'tlattend','tlignore','atttrials','igntrials')
    %%%%%
  end
  
  
  %%%%% get tl for pitch binaural
  % determine which stimuli were lo/hi
  switch info.(currsubj).leftrightelec
    case 'right'
      lotrials = attBLlostd;
      hitrials = attBLhistd;
      
    case 'left'
      lotrials = attBRlostd;
      hitrials = attBRhistd;
      
    otherwise
      error('woopsie')
  end
  % get timelocked average
  cfgtl.trials           = lotrials;
  tllopitch = ft_timelockanalysis(cfgtl,filtdata);
  cfgtl.trials           = hitrials;
  tlhipitch = ft_timelockanalysis(cfgtl,filtdata);
  % save shit
  fn = [savepath currsubj '_' 'HFA100-10-200_pitchbin' '.mat'];
  save(fn,'tllopitch','tlhipitch','lotrials','hitrials')
  %%%%%
  
  
  
  
  %%%%% get tl for pitch max trials
  % determine which stimuli were lo/hi
  % determine which stimuli were lo/hi while ensuring equal trial numbers for lo/hi
  lotrials = [];
  hitrials = [];
  % attRR
  dum = min([numel(attRRlostd) numel(attRRhistd)]);
  lotrials = [lotrials; attRRlostd(1:dum)];
  hitrials = [hitrials; attRRhistd(1:dum)];
  % attRL
  dum = min([numel(attRLlostd) numel(attRLhistd)]);
  lotrials = [lotrials; attRLlostd(1:dum)];
  hitrials = [hitrials; attRLhistd(1:dum)];
  % attLR
  dum = min([numel(attLRlostd) numel(attLRhistd)]);
  lotrials = [lotrials; attLRlostd(1:dum)];
  hitrials = [hitrials; attLRhistd(1:dum)];
  % attLL
  dum = min([numel(attLLlostd) numel(attLLhistd)]);
  lotrials = [lotrials; attLLlostd(1:dum)];
  hitrials = [hitrials; attLLhistd(1:dum)];
  % attBR
  dum = min([numel(attBRlostd) numel(attBRhistd)]);
  lotrials = [lotrials; attBRlostd(1:dum)];
  hitrials = [hitrials; attBRhistd(1:dum)];
  % attBL
  dum = min([numel(attBLlostd) numel(attBLhistd)]);
  lotrials = [lotrials; attBLlostd(1:dum)];
  hitrials = [hitrials; attBLhistd(1:dum)];
  
  % get timelocked average
  cfgtl.trials           = lotrials;
  tllopitch = ft_timelockanalysis(cfgtl,filtdata);
  cfgtl.trials           = hitrials;
  tlhipitch = ft_timelockanalysis(cfgtl,filtdata);
  % save shit
  fn = [savepath currsubj '_' 'HFA100-10-200_pitchmaxtrials' '.mat'];
  save(fn,'tllopitch','tlhipitch','lotrials','hitrials')
  %%%%%
  
  


  
  
  
  
  
end






























function playground




% get details
info = rmr_abcaudecog_info;

% set stuff
hfapath   = [info.savepath 'hfa/'];

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  
  % get layout
  load([info.savepath currsubj '_fieldtrip_layout.mat']);
  
  
  % plot tl for attention
  for iset = 1:numel(info.(currsubj).lohipitch)
    fn = [hfapath currsubj '_' 'HFA100-10-200_attign' '_' info.(currsubj).lohipitch{iset} '-pitch' '.mat'];
    load(fn)
    
    % make time axis identical
    ind = 1:min(numel(tlattend.time),numel(tlignore.time));
    tlattend.avg    = tlattend.avg(:,ind);
    tlattend.var    = tlattend.var(:,ind);
    tlattend.time   = tlattend.time(ind);
    tlattend.dof    = tlattend.dof(:,ind);
    tlattend.trial  = tlattend.trial(:,:,ind);
    tlignore.avg    = tlignore.avg(:,ind);
    tlignore.var    = tlignore.var(:,ind);
    tlignore.time   = tlignore.time(ind);
    tlignore.dof    = tlignore.dof(:,ind);
    tlignore.trial  = tlignore.trial(:,:,ind);
      
    % plot 
    figure('name',[currsubj ' att' '_' info.(currsubj).lohipitch{iset}],'numbertitle','off')
    cfg = [];
    cfg.layout      = lay;
    cfg.showoutline = 'yes';
    ft_multiplotER(cfg,tlignore,tlattend)
    title(['ign(b) att(r)       ntrial ' num2str(size(tlignore.trial,1)) '|' num2str(size(tlattend.trial,1))])
    % calc t and plot
    cfg = [];
    cfg.method    = 'analytic';
    cfg.statistic = 'ft_statfun_indepsamplesT';
    cfg.design    = [ones(1,size(tlattend.trial,1)) ones(1,size(tlignore.trial,1))*2];
    cfg.ivar      = 1;
    tltstat = ft_timelockstatistics(cfg,tlignore,tlattend);
    title(['ign(b) att(r)       ntrial ' num2str(size(tlignore.trial,1)) '|' num2str(size(tlattend.trial,1))])
    figure
    cfg = [];
    cfg.layout      = lay;
    cfg.showoutline = 'yes';
    cfg.parameter   = 'stat';
    cfg.vlim        = [-5 5];
    ft_multiplotER(cfg,tltstat)
    title(['ign(b) att(r)       ntrial ' num2str(size(tlignore.trial,1)) '|' num2str(size(tlattend.trial,1))])
  end
  
  
  % plot tl for pitch
  for iplot = 1:2
    if iplot==1
      fn = [hfapath currsubj '_' 'HFA100-10-200_pitchbin' '.mat'];
    elseif iplot==2
      fn = [hfapath currsubj '_' 'HFA100-10-200_pitchmaxtrials' '.mat'];
    end
    load(fn)
    
    % make time axis identical
    ind = 1:min(numel(tllopitch.time),numel(tlhipitch.time));
    tllopitch.avg    = tllopitch.avg(:,ind);
    tllopitch.var    = tllopitch.var(:,ind);
    tllopitch.time   = tllopitch.time(ind);
    tllopitch.dof    = tllopitch.dof(:,ind);
    tllopitch.trial  = tllopitch.trial(:,:,ind);
    tlhipitch.avg    = tlhipitch.avg(:,ind);
    tlhipitch.var    = tlhipitch.var(:,ind);
    tlhipitch.time   = tlhipitch.time(ind);
    tlhipitch.dof    = tlhipitch.dof(:,ind);
    tlhipitch.trial  = tlhipitch.trial(:,:,ind);
    
    % plot
    figure
    cfg = [];
    cfg.layout      = lay;
    cfg.showoutline = 'yes';
    ft_multiplotER(cfg,tllopitch,tlhipitch)
    title(['lo(b) hi(r)       ntrial ' num2str(size(tllopitch.trial,1)) '|' num2str(size(tlhipitch.trial,1))])
    % calc t and plot
    cfg = [];
    cfg.method    = 'analytic';
    cfg.statistic = 'ft_statfun_indepsamplesT';
    cfg.design    = [ones(1,size(tllopitch.trial,1)) ones(1,size(tlhipitch.trial,1))*2];
    cfg.ivar      = 1;
    tltstat = ft_timelockstatistics(cfg,tllopitch,tlhipitch);
    title(['lo(b) hi(r)       ntrial ' num2str(size(tllopitch.trial,1)) '|' num2str(size(tlhipitch.trial,1))])
    figure
    cfg = [];
    cfg.layout      = lay;
    cfg.showoutline = 'yes';
    cfg.parameter   = 'stat';
    ft_multiplotER(cfg,tltstat)
    title(['lo(b) hi(r)       ntrial ' num2str(size(tllopitch.trial,1)) '|' num2str(size(tlhipitch.trial,1))])
  end
end
























% get details
info = rmr_abcaudecog_info;

% set stuff
hfapath   = [info.savepath 'hfa/'];

% sub select, get rid of the JH's
info.subj = {'GP15','GP22','GP28','GP35','ST1','ST6','ST8'};

% sub select for hfa differences between conditions of some sort
%info.subj = {'GP15','GP22','GP28','GP35'};

% set glboal axis limits
vlim = [-.5 .5];

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  % get layout
  load([info.savepath currsubj '_fieldtrip_layout.mat']);
  
  
  %%% PLOT ATTENTION
  for iset = 1:numel(info.(currsubj).lohipitch)
    fn = [hfapath currsubj '_' 'HFA100-10-200_attign' '_' info.(currsubj).lohipitch{iset} '-pitch' '.mat'];
    load(fn)
    
    % make time axis identical
    ind = 1:min(numel(tlattend.time),numel(tlignore.time));
    tlattend.avg    = tlattend.avg(:,ind);
    tlattend.var    = tlattend.var(:,ind);
    tlattend.time   = tlattend.time(ind);
    tlattend.dof    = tlattend.dof(:,ind);
    tlattend.trial  = tlattend.trial(:,:,ind);
    tlignore.avg    = tlignore.avg(:,ind);
    tlignore.var    = tlignore.var(:,ind);
    tlignore.time   = tlignore.time(ind);
    tlignore.dof    = tlignore.dof(:,ind);
    tlignore.trial  = tlignore.trial(:,:,ind);
    
    % get patch coordinates
    avg1 = tlignore.avg;
    sem1 = sqrt(tlignore.var)./(sqrt(size(tlignore.trial,1)));
    avg2 = tlattend.avg;
    sem2 = sqrt(tlattend.var)./(sqrt(size(tlattend.trial,1)));
    time = tlattend.time;
    %vlim = [round(min(min(min([avg1-sem1 avg2-sem2])))*10)./10 round(max(max(max([avg1+sem1 avg2+sem2])))*10)./10];
    
    % remove NaNs
    nanind = isnan(sum(avg1+avg2+sem1+sem2,1));
    avg1 = avg1(:,~nanind);
    avg2 = avg2(:,~nanind);
    sem1 = sem1(:,~nanind);
    sem2 = sem2(:,~nanind);
    time = time(:,~nanind);
    
    % plot the patch
    figure('name',[currsubj ' att' '_' info.(currsubj).lohipitch{iset}],'numbertitle','off')
    hold on
    % draw outlines
    for ioutline = 1:numel(lay.outline)
      xline = lay.outline{ioutline}(:,1);
      yline = lay.outline{ioutline}(:,2);
      line(xline,yline,'linewidth',1,'linestyle','--')
    end
    % draw data
    label = tlattend.label;
    nchan = numel(label);
    for ichan = 1:nchan
      % get local axis
      layind = strcmp(lay.label,label{ichan});
      hpos   = lay.pos(layind,1);
      vpos   = lay.pos(layind,2);
      width  = lay.width(layind);
      height = lay.height(layind);
      
      % plot patch
      x = [time time(end:-1:1)];
      y = [avg1(ichan,:)-sem1(ichan,:) avg1(ichan,end:-1:1)+sem1(ichan,end:-1:1)];
      ft_plot_patch(x,y,'hpos',hpos,'vpos',vpos,'height',height,'width',width,'axis','yes','vlim',vlim,'facecolor','b','edgecolor','none','facealpha',.5)
      y = [avg2(ichan,:)-sem2(ichan,:) avg2(ichan,end:-1:1)+sem2(ichan,end:-1:1)];
      ft_plot_patch(x,y,'hpos',hpos,'vpos',vpos,'height',height,'width',width,'axis','yes','vlim',vlim,'facecolor','r','edgecolor','none','facealpha',.5)
    end
    axis off
    title(['ign(b) att(r)   ylim ' num2str(vlim(1)) ' - ' num2str(vlim(2)) '    ntrial ' num2str(size(tlignore.trial,1)) '|' num2str(size(tlattend.trial,1))])
  end
  
  
  
  %%% PLOT PITCH
  for iplot = 1:2
    if iplot==1
      fn = [hfapath currsubj '_' 'HFA100-10-200_pitchbin' '.mat'];
    elseif iplot==2
      fn = [hfapath currsubj '_' 'HFA100-10-200_pitchmaxtrials' '.mat'];
    end
    load(fn)
    
    % make time axis identical
    ind = 1:min(numel(tllopitch.time),numel(tlhipitch.time));
    tllopitch.avg    = tllopitch.avg(:,ind);
    tllopitch.var    = tllopitch.var(:,ind);
    tllopitch.time   = tllopitch.time(ind);
    tllopitch.dof    = tllopitch.dof(:,ind);
    tllopitch.trial  = tllopitch.trial(:,:,ind);
    tlhipitch.avg    = tlhipitch.avg(:,ind);
    tlhipitch.var    = tlhipitch.var(:,ind);
    tlhipitch.time   = tlhipitch.time(ind);
    tlhipitch.dof    = tlhipitch.dof(:,ind);
    tlhipitch.trial  = tlhipitch.trial(:,:,ind);
    
    % get patch coordinates
    avg1 = tllopitch.avg;
    sem1 = sqrt(tllopitch.var)./(sqrt(size(tllopitch.trial,1)));
    avg2 = tlhipitch.avg;
    sem2 = sqrt(tlhipitch.var)./(sqrt(size(tlhipitch.trial,1)));
    time = tllopitch.time;
    %vlim = [round(min(min(min([avg1-sem1 avg2-sem2])))*10)./10 round(max(max(max([avg1+sem1 avg2+sem2])))*10)./10];
    vlim = [-.5 .5];
    % remove NaNs
    nanind = isnan(sum(avg1+avg2+sem1+sem2,1));
    avg1 = avg1(:,~nanind);
    avg2 = avg2(:,~nanind);
    sem1 = sem1(:,~nanind);
    sem2 = sem2(:,~nanind);
    time = time(:,~nanind);
    
    
    % plot the patch
    figure('name',[currsubj ' pitch'],'numbertitle','off')
    hold on
    % draw outlines
    for ioutline = 1:numel(lay.outline)
      xline = lay.outline{ioutline}(:,1);
      yline = lay.outline{ioutline}(:,2);
      line(xline,yline,'linewidth',1,'linestyle','--')
    end
    % draw data
    label = tlattend.label;
    nchan = numel(label);
    for ichan = 1:nchan
      % get local axis
      layind = strcmp(lay.label,label{ichan});
      hpos   = lay.pos(layind,1);
      vpos   = lay.pos(layind,2);
      width  = lay.width(layind);
      height = lay.height(layind);
      
      % plot patch
      x = [time time(end:-1:1)];
      y = [avg1(ichan,:)-sem1(ichan,:) avg1(ichan,end:-1:1)+sem1(ichan,end:-1:1)];
      ft_plot_patch(x,y,'hpos',hpos,'vpos',vpos,'height',height,'width',width,'axis','yes','vlim',vlim,'facecolor','b','edgecolor','none','facealpha',.5)
      y = [avg2(ichan,:)-sem2(ichan,:) avg2(ichan,end:-1:1)+sem2(ichan,end:-1:1)];
      ft_plot_patch(x,y,'hpos',hpos,'vpos',vpos,'height',height,'width',width,'axis','yes','vlim',vlim,'facecolor','r','edgecolor','none','facealpha',.5)
    end
    axis off
    title(['lo(b) hi(r)   ylim ' num2str(vlim(1)) ' - ' num2str(vlim(2)) '    ntrial ' num2str(size(tllopitch.trial,1)) '|' num2str(size(tlhipitch.trial,1))])
  end
end










