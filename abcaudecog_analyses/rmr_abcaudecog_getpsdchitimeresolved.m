function rmr_abcaudecog_getpsdchitimeresolved

warning('experimental!!!!')

% get details
info = rmr_abcaudecog_info;

% set stuff
psdpath   = [info.savepath 'psd/'];
condlabel = {'att','ign'};

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1 %   :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  % skip subj if done already
  fn{1} = [psdpath currsubj '_' condlabel{1} '_75msoffschi_timeres_80-500.mat'];
  fn{2} = [psdpath currsubj '_' condlabel{2} '_75msoffschi_timeres_80-500.mat'];
  if exist(fn{1},'file') && exist(fn{2},'file')
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
  
  % get chi and offset
  condtrials = {atttrials,igntrials};
  for icond = 1:2
    
    % get freq with time
    cfg = [];
    cfg.channel    = 'all';
    cfg.trials     = condtrials{icond};
    cfg.keeptrials = 'yes';
    cfg.pad        = 1; % pad out to 1 second;
    cfg.method     = 'mtmconvol';
    cfg.output     = 'pow';
    cfg.toi        = toi;
    cfg.foi        = 80:500;
    cfg.taper      = 'hanning';
    cfg.t_ftimwin = ones(numel(cfg.foi),1) .* .075;
    freq = ft_freqanalysis(cfg,data);
    
    % log the freq and pow
    freq.powspctrm = log10(freq.powspctrm);
    freq.freq = log10(freq.freq);
    
    % extract chi (WITHOUT MINUS SIGN!)
    % fit ax = b using robust fit
    nchan = numel(freq.label);
    nrpt  = size(freq.powspctrm,1);
    ntime = size(freq.powspctrm,4);
    offschi = NaN(nrpt,nchan,2,ntime);
    powat80 = NaN(nrpt,nchan,ntime);
    for itrial = 1:size(freq.powspctrm,1)
      disp(['fitting offset/chi for trial #' num2str(itrial)])
      for ichan = 1:nchan
        for itime = 1:ntime
          % att
          y = squeeze(freq.powspctrm(itrial,ichan,:,itime));
          if any(isnan(y))
            continue
          end
          x = freq.freq;
          powat80(itrial,ichan,itime) = y(1); % save pow at 80Hz
          offschi(itrial,ichan,:,itime) = rmr_robustfit(x,y);
        end
      end
    end
    offschi(:,:,2,:) = offschi(:,:,2,:) .* -1;
       
    % save and clear
    save(fn{icond} ,'offschi','freq','powat80')
    clear freq offschi
  end
end



function playground


















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
  fn{1} = [psdpath currsubj '_' condlabel{1} '_75msoffschi_timeres_80-500.mat'];
  fn{2} = [psdpath currsubj '_' condlabel{2} '_75msoffschi_timeres_80-500.mat'];

  
  % plot
  for icond = 1:2
    load(fn{icond});
    
    % calc mean and std of chi
    meanchi = squeeze(mean(offschi(:,:,2,:),1));
    semchi  = squeeze(std(offschi(:,:,2,:),1));%squeeze(std(offschi(:,:,2,:),1))./sqrt(size(offschi,1));
    vlim = [0 round(max(max(max(meanchi)))*10)./10];

    % plot local pac
    figure('name',[currsubj ' ' condlabel{icond}])
    hold on
    % draw outlines
    for ioutline = 1:numel(lay.outline)
      xline = lay.outline{ioutline}(:,1);
      yline = lay.outline{ioutline}(:,2);
      line(xline,yline,'linewidth',1,'linestyle','--')
    end
    % draw data
    label = freq.label;
    nchan = numel(label);
    for ichan = 1:nchan
      % get local axis
      layind = strcmp(lay.label,label{ichan});
      hpos   = lay.pos(layind,1);
      vpos   = lay.pos(layind,2);
      width  = lay.width(layind);
      height = lay.height(layind);
      
      % plot
      x = freqwphase.time_used;
      y = cat(1,meanchi(ichan,:),meanchi(ichan,:)-semchi(ichan,:),meanchi(ichan,:)+semchi(ichan,:));     
      ft_plot_vector(x,y,'hpos',hpos,'vpos',vpos,'height',height,'width',width,'color','b','axis','yes','vlim',vlim)
      hline = findobj(gca,'type','line');
      set(hline(3),'color',[.5 .5 1])
      set(hline(4),'color',[.5 .5 1])
    end
    axis off
    title(['ylim ' num2str(vlim(1)) ' - ' num2str(vlim(2))])
  end
  
  
  
end








