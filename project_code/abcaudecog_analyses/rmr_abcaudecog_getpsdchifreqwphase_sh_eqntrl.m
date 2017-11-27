function rmr_abcaudecog_getpsdchifreqwphase_sh_eqntrl


% get details
info = rmr_abcaudecog_info;

% set stuff
psdpath   = [info.savepath 'psd/'];

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1  %  :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  % skip subj if done already
  fnfreqwphase = [psdpath currsubj '_' 'att-ign' '_freqwphase_sh_eqntrl_3cyc75msoffschi80-500.mat'];
  if exist(fnfreqwphase,'file')
    continue
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
    cfg.bsfilter   = 'yes';
    cfg.bsfilttype = 'but';
    cfg.bsfiltord  = 2;
    cfg.bsfreq     = info.(currsubj).bsfreq;
    cfg.lpfilter   = 'yes';
    cfg.lpfreq     = 540;
    cfg.lpfilttype = 'fir';
    cfg.lpfiltord  = 1000;
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
  else
    data = datacmb{1};
  end
  
  
  % set trial bool
  attRRstd = data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4; % attend right, right stim, standard,  correctly rejected
  attRLstd = data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4; % attend right, left  stim, standard,  correctly rejected
  attLRstd = data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 2 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4; % attend left,  right stim, standard,  correctly rejected
  attLLstd = data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 1 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4; % attend left,  left stim,  standard,  correctly rejected
  attLBstd = data.trialinfo(:,1) == 1 & data.trialinfo(:,6) == 3 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4; % attend left,  binau stim, standard,  correctly rejected
  attRBstd = data.trialinfo(:,1) == 2 & data.trialinfo(:,6) == 3 & data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4; % attend right, binau stim, standard,  correctly rejected
  lopitch = data.trialinfo(:,7) == 1;
  hipitch = data.trialinfo(:,7) == 2;
  
  % determine which stimuli were (un)attended
  switch info.(currsubj).leftrightelec
    case 'right'
      atttrials = find(attRRstd);
      igntrials = find(attLRstd);
    case 'left'
      atttrials = find(attLLstd);
      igntrials = find(attRLstd);
    otherwise
      error('woopsie')
  end
  % set toi sampling
  if data.fsample>3000
    [dum maxind] = max(cellfun(@numel,data.time));
    toi = data.time{maxind}(1:8:end);
  else
    [dum maxind] = max(cellfun(@numel,data.time));
    toi = data.time{maxind}(1:2:end);
  end
  
  % downselect trials if one of the conditions is bigger, ALSO USED for offschi below
  disp(['randomly equalizing ntrials of ' num2str(numel(atttrials)) ' and ' num2str(numel(igntrials))])
  trlind = [];
  if numel(atttrials) < numel(igntrials)
    trlind{2} = randperm(numel(igntrials));
    trlind{2} = trlind{2}(1:numel(atttrials));
    igntrials = igntrials(trlind{2});
    trlind{1} = 1:numel(trlind{2});
  elseif numel(igntrials) < numel(atttrials)
    trlind{1} = randperm(numel(atttrials));
    trlind{1} = trlind{1}(1:numel(igntrials));
    atttrials = atttrials(trlind{1});
    trlind{2} = 1:numel(trlind{1});
  end
  
  
  % get wPLFs
  condtrials = {atttrials,igntrials};
  condlabel  = {'att','ign'};
  freqwphaseattign    = [];
  freqwphaseattignsh1 = [];
  freqwphaseattignsh2 = [];
  for icond = 1:2
    % fourfreq
    cfg = [];
    cfg.channel    = 'all';
    cfg.trials     = condtrials{icond};
    cfg.keeptrials = 'yes';
    cfg.pad        = 1; % pad out to 1 second;
    cfg.method     = 'mtmconvol';
    cfg.output     = 'fourier';
    cfg.toi        = toi;
    cfg.foi        = 8:30;
    cfg.t_ftimwin  = 3./cfg.foi;
    cfg.taper      = 'hanning';
    freqphas = ft_freqanalysis(cfg,data);
    
    % chi and offset if not on disk
    fnoffschi = [psdpath currsubj '_' condlabel{icond} '_75msoffschi_timeres_80-500.mat'];
    load(fnoffschi)
    offschi = offschi(trlind{icond},:,:,:);
    
    
    % create freqdata
    freqdata = rmfield(freqphas,{'fourierspctrm','freq'});
    freqdata.phasfreqfour = freqphas.fourierspctrm;
    freqdata.phasfreq     = freqphas.freq;
    freqdata.ampfreqenv   = offschi;
    freqdata.ampfreq      = [1 2]; % include both to be sure, in case of silly n=1 errors for some of the dimensions in roe_calc...
    freqdata.dimord       = 'rpt_chan_freq_time';
    clear freqfour
    
    % calc freqwphase for all trials
    cfg = [];
    cfg.phasparam = 'phasfreqfour';
    cfg.ampparam  = 'ampfreqenv';
    freqwphaseattign{icond} = roe_calcwphaselocking(cfg,freqdata);
    
    % calc freqwphase for odd even split
    % odd
    freqdatash1 = freqdata;
    freqdatash1.phasfreqfour = freqdatash1.phasfreqfour(1:2:end,:,:,:);
    freqdatash1.ampfreqenv   = freqdatash1.ampfreqenv(1:2:end,:,:,:);
    cfg = [];
    cfg.phasparam = 'phasfreqfour';
    cfg.ampparam  = 'ampfreqenv';
    freqwphaseattignsh1{icond} = roe_calcwphaselocking(cfg,freqdatash1);
    % even
    freqdatash2 = freqdata;
    freqdatash2.phasfreqfour = freqdatash2.phasfreqfour(2:2:end,:,:,:);
    freqdatash2.ampfreqenv   = freqdatash2.ampfreqenv(2:2:end,:,:,:);
    cfg = [];
    cfg.phasparam = 'phasfreqfour';
    cfg.ampparam  = 'ampfreqenv';
    freqwphaseattignsh2{icond} = roe_calcwphaselocking(cfg,freqdatash2);
  end % icond
  
  
  
  % combine wphaselocks of all of the above
  % all trials
  freqwphase = rmfield(freqwphaseattign{1},'ampfreq');
  freqwphase.wphaselockspctrm = cat(4,squeeze(freqwphaseattign{1}.wphaselockspctrm(:,:,2,:)),squeeze(freqwphaseattign{2}.wphaselockspctrm(:,:,2,:)));
  freqwphase.dimord   = 'ampchan_phaschan_phasfreq_cond';
  freqwphase.label    = freqwphase.label([1 2 4]);
  freqwphase.label{4} = {'att','ign'};
  freqwphase.trial_used = cat(1,freqwphaseattign{1}.trials_used,freqwphaseattign{2}.trials_used);
  freqwphase.time_used  = cat(1,freqwphaseattign{1}.time_used,freqwphaseattign{2}.time_used);
  freqwphase.avgtimecnt = permute(cat(3,freqwphaseattign{1}.avgtimecnt,freqwphaseattign{2}.avgtimecnt),[3 1 2]);
  
  % sh1
  freqwphasesh1 = rmfield(freqwphaseattignsh1{1},'ampfreq');
  freqwphasesh1.wphaselockspctrm = cat(4,squeeze(freqwphaseattignsh1{1}.wphaselockspctrm(:,:,2,:)),squeeze(freqwphaseattignsh1{2}.wphaselockspctrm(:,:,2,:)));
  freqwphasesh1.dimord   = 'ampchan_phaschan_phasfreq_cond';
  freqwphasesh1.label    = freqwphasesh1.label([1 2 4]);
  freqwphasesh1.label{4} = {'att','ign'};
  freqwphasesh1.trial_used = cat(1,freqwphaseattignsh1{1}.trials_used,freqwphaseattignsh1{2}.trials_used);
  freqwphasesh1.time_used  = cat(1,freqwphaseattignsh1{1}.time_used,freqwphaseattignsh1{2}.time_used);
  freqwphasesh1.avgtimecnt = permute(cat(3,freqwphaseattignsh1{1}.avgtimecnt,freqwphaseattignsh1{2}.avgtimecnt),[3 1 2]);
  
  % sh2
  freqwphasesh2 = rmfield(freqwphaseattignsh2{1},'ampfreq');
  freqwphasesh2.wphaselockspctrm = cat(4,squeeze(freqwphaseattignsh2{1}.wphaselockspctrm(:,:,2,:)),squeeze(freqwphaseattignsh2{2}.wphaselockspctrm(:,:,2,:)));
  freqwphasesh2.dimord   = 'ampchan_phaschan_phasfreq_cond';
  freqwphasesh2.label    = freqwphasesh2.label([1 2 4]);
  freqwphasesh2.label{4} = {'att','ign'};
  freqwphasesh2.trial_used = cat(1,freqwphaseattignsh2{1}.trials_used,freqwphaseattignsh2{2}.trials_used);
  freqwphasesh2.time_used  = cat(1,freqwphaseattignsh2{1}.time_used,freqwphaseattignsh2{2}.time_used);
  freqwphasesh2.avgtimecnt = permute(cat(3,freqwphaseattignsh2{1}.avgtimecnt,freqwphaseattignsh2{2}.avgtimecnt),[3 1 2]);
  
  
  % save and clear
  save(fnfreqwphase,'freqwphase','freqwphasesh1','freqwphasesh2','offschi','atttrials','igntrials')
  clear freqdata freqwphase freqwphasesh1 freqwphasesh2 freqwphaseattign freqwphaseattignsh1 freqwphaseattignsh2
end























function playground










% get details
info = rmr_abcaudecog_info;

% set stuff
psdpath   = [info.savepath 'psd/'];

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1  %  :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  fnfreqwphase = [psdpath currsubj '_' 'att-ign' '_freqwphase_sh_eqntrl_3cyc75msoffschi80-500.mat'];
  load(fnfreqwphase)
  
  % prep
  load([info.savepath currsubj '_fieldtrip_layout.mat']);
  % remove channels from lay
  selchan = match_str(lay.label,freqwphase.label_old);
  lay.pos    = lay.pos(selchan,:);
  lay.label  = lay.label(selchan);
  lay.width  = lay.width(selchan);
  lay.height = lay.height(selchan);
  
  
  % browse, but trick the browse code first
  freqwphase.wphaselockspctrm = permute(freqwphase.wphaselockspctrm,[1 2 4 3]);
  freqwphase.ampfreq = [1 2];
  cfg = [];
  cfg.layout     = lay;
  cfg.windowname = [currsubj];
  roe_browsefreqwphase(cfg,freqwphase)
  
end % isubj




















% get details
info = rmr_abcaudecog_info;

% set paths
psdpath   = [info.savepath 'psd/'];
condlabel = {'att','ign'};

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1   :nsubj
  
  % set
  currsubj = info.subj{isubj};
  load([info.savepath currsubj '_fieldtrip_layout.mat']);
  
  fnfreqwphase = [psdpath currsubj '_' 'att-ign' '_freqwphase_sh_eqntrl_3cyc75msoffschi80-500.mat'];
  load(fnfreqwphase)
  
  % select WITHIN CHAN PAC ONLY
  nchan = numel(freqwphase.label_old);
  newwplfsize = size(freqwphase.wphaselockspctrm);
  newwplfsize = newwplfsize([2 3 4]);
  freqwphase.newwphaselockspctrm = NaN(newwplfsize);
  for ichan = 1:nchan
    freqwphase.newwphaselockspctrm(ichan,:,:)    = squeeze(freqwphase.wphaselockspctrm(ichan,ichan,:,:));
  end
  freqwphase.wphaselockspctrm    = freqwphase.newwphaselockspctrm;
  freqwphase    = rmfield(freqwphase,'newwphaselockspctrm');
  
  
  % set clim based on both conditions
  wplf = abs(freqwphase.wphaselockspctrm);
  maxwplf = round(max(max(max(abs(freqwphase.wphaselockspctrm))))*100)/100;
  nchan = numel(freqwphase.label_old);
  maxwplf = .1;
  % plot local pac
  figure('name',[currsubj])
  hold on
  % draw outlines
  for ioutline = 1:numel(lay.outline)
    xline = lay.outline{ioutline}(:,1);
    yline = lay.outline{ioutline}(:,2);
    line(xline,yline,'linewidth',1,'linestyle','--')
  end
  % draw data
  label = freqwphase.label_old;
  nchan = numel(freqwphase.label_old);
  for ichan = 1:nchan
    % get local axis
    layind = strcmp(lay.label,label{ichan});
    hpos   = lay.pos(layind,1);
    vpos   = lay.pos(layind,2);
    width  = lay.width(layind);
    height = lay.height(layind);
    
    % plot
    c = abs(squeeze(freqwphase.wphaselockspctrm(ichan,:,:))).';
    % generated fitted line and plot
    y = freqwphase.phasfreq;
    x = [1 2];
    clim = [0  maxwplf];
    ft_plot_matrix(x,y,c,'hpos',hpos,'vpos',vpos,'height',height,'width',width,'color','brbr','clim',clim)
  end
  axis off
  title(['clim ' num2str(clim(1)) ' - ' num2str(clim(2))])
  
   
end



























