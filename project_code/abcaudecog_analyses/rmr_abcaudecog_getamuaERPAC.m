function rmr_abcaudecog_getamuaERPAC



% get details
info = rmr_abcaudecog_info;

% set stuff
savepath   = [info.savepath 'amua/'];

% sub select for amua differences between conditions of some sort
info.subj = {'GP15','GP22','GP28','GP35'};

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  % check for already done
  if numel(info.(currsubj).lohipitch)==2
    fn1 = [savepath currsubj '_' 'AMUA100-10-200_ERPAC_attign' '_' info.(currsubj).lohipitch{1} '-pitch' '.mat'];
    fn2 = [savepath currsubj '_' 'AMUA100-10-200_ERPAC_attign' '_' info.(currsubj).lohipitch{2} '-pitch' '.mat'];
    fn3 = [savepath currsubj '_' 'AMUA100-10-200_ERPAC_pitch' '.mat'];
    if all([exist(fn1,'file') exist(fn2,'file') exist(fn3,'file')])
      continue
    end
  else
    fn1 = [savepath currsubj '_' 'AMUA100-10-200_ERPAC_attign' '_' info.(currsubj).lohipitch{1} '-pitch' '.mat'];
    fn2 = [savepath currsubj '_' 'AMUA100-10-200_ERPAC_pitch' '.mat'];
    if all([exist(fn1,'file') exist(fn2,'file')])
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
    cfg.bsfilter   = 'yes';
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
    % add 300ms baseline
    trl(:,1) = trl(:,1) - round(0.3 .* data.fsample);
    trl(:,3) = -round(0.3 .* data.fsample);
    %%%%%%%%%%%
    % ADD 200ms POST TRIAL
    trl(:,2) = trl(:,2) + round(0.2 .* data.fsample);
    %%%%%%%%%%%
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
  
  % set toi sampling
  [dum maxind] = max(cellfun(@numel,data.time));
  toi = data.time{maxind}(round(data.fsample.*.1):1:end-round(.3.*data.fsample)); % remove the extra ms from toi, and throw away another 100ms to be sure 
  if data.fsample>3000
    [dum maxind] = max(cellfun(@numel,data.time));
    toi = toi(1:4:end);
  else
    % no need
    %         [dum maxind] = max(cellfun(@numel,data.time));
    %         toi = toi(1:2:end);
  end
  
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
    disp(['getting AMUA ampspctrm of ' num2str(freqoi(ifreq)) 'Hz'])
    % get powspctrm
    cfg = [];
    cfg.feedback   = 'none';
    cfg.channel    = 'all';
    cfg.keeptrials = 'yes';
    cfg.pad        = 2; % pad out to 1 second;
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
  %   % remove NaNs from filtdata (raw data cannot contain NaNs)
  %   for itrial = 1:ntrial
  %     nanindfiltdata = isnan(sum(filtdata.trial{itrial},1));
  %     filtdata.trial{itrial} = filtdata.trial{itrial}(:,~nanindfiltdata);
  %     filtdata.time{itrial}  = filtdata.time{itrial}(:,~nanindfiltdata);
  %   end
  clear freqdata ampenvspctrm tmp
  
 
  %%%%% get PAC for attention
  for iset = 1:numel(info.(currsubj).lohipitch)
    % determine which stimuli were (un)attended
    switch info.(currsubj).leftrightelec
      case 'right'
        if      strcmp(info.(currsubj).lohipitch{iset},'low')
          atttrials = find(attRRlostd);
          igntrials = find(attLRlostd);
        elseif  strcmp(info.(currsubj).lohipitch{iset},'high')
          atttrials = find(attRRhistd);
          igntrials = find(attLRhistd);
        else
          error('woopsie')
        end
        
      case 'left'
        if     strcmp(info.(currsubj).lohipitch{iset},'low')
          atttrials = find(attLLlostd);
          igntrials = find(attRLlostd);
        elseif strcmp(info.(currsubj).lohipitch{iset},'high')
          atttrials = find(attLLhistd);
          igntrials = find(attRLhistd);
        else
          error('woopsie')
        end
        
      otherwise
        error('woopsie')
    end
        
    % get fourfreq
    cfg = [];
    cfg.channel    = 'all';
    cfg.trials     = [atttrials' igntrials'];
    cfg.keeptrials = 'yes';
    cfg.pad        = 2; % pad out to 1 second;
    cfg.method     = 'mtmconvol';
    cfg.output     = 'fourier';
    cfg.toi        = toi; 
    cfg.foi        = 8:30;
    cfg.t_ftimwin  = 3./cfg.foi;
    cfg.taper      = 'hanning';
    freqphas = ft_freqanalysis(cfg,data);

    % get amua array
    amua = filtdata.trial([atttrials' igntrials']);
    amua = permute(cat(3,amua{:}),[3 1 2]);
    amua = cat(4,amua,amua);
    amua = permute(amua,[1 2 4 3]);
       
    % create freqdata
    freqdata = rmfield(freqphas,{'fourierspctrm','freq'});
    freqdata.phasfreqfour = freqphas.fourierspctrm;
    freqdata.phasfreq     = freqphas.freq;
    freqdata.ampfreqenv   = amua;
    freqdata.ampfreq      = [1 2]; % include both to be sure, in case of silly n=1 errors for some of the dimensions in roe_calc...
    freqdata.dimord       = 'rpt_chan_freq_time';
    clear freqphas
    
    % pre remove NaNs
    nanampenv    = isnan(squeeze(mean(mean(mean(freqdata.ampfreqenv,1),2),3)));
    nanfourierts = isnan(squeeze(mean(mean(mean(freqdata.phasfreqfour,1),2),3)));
    nonnantimepoints = logical(~nanampenv .* ~nanfourierts);
    freqdata.phasfreqfour = freqdata.phasfreqfour(:,:,:,nonnantimepoints);
    freqdata.ampfreqenv   = freqdata.ampfreqenv(:,:,:,nonnantimepoints);  
    display(['new # time-points = ' num2str(sum(nonnantimepoints)) ' (removed ' num2str(numel(toi) - sum(nonnantimepoints)) ' timepoints)'])
    % save time axis to put back later
    timeinfo         = freqdata.time(nonnantimepoints);
    
    % split it
    % for att
    freqatt = freqdata;
    freqatt.phasfreqfour = permute(freqatt.phasfreqfour(1:numel(atttrials),:,:,:), [4 2 3 1]);
    freqatt.ampfreqenv   = permute(freqatt.ampfreqenv(1:numel(atttrials),:,:,:), [4 2 3 1]);
    freqatt.time         = atttrials;
    % for ign
    freqign = freqdata;
    freqign.phasfreqfour = permute(freqign.phasfreqfour(numel(atttrials)+1:end,:,:,:), [4 2 3 1]);
    freqign.ampfreqenv   = permute(freqign.ampfreqenv(numel(atttrials)+1:end,:,:,:), [4 2 3 1]);
    freqign.time         = igntrials;
    clear freqdata

    % calc freqwphase and insert timeinfo
    cfg = [];
    cfg.handlenans    = 'conservative';
    cfg.normalization = 'withintrials';
    cfg.phasparam     = 'phasfreqfour';
    cfg.ampparam      = 'ampfreqenv';
    cfg.keeptrials    = 'yes';
    % for att
    attamuapac = roe_calcwphaselocking(cfg,freqatt);
    attamuapac.timeinfo = timeinfo;
    % for ign
    ignamuapac = roe_calcwphaselocking(cfg,freqign);
    ignamuapac.timeinfo = timeinfo;
      
    % save shit
    fn = [savepath currsubj '_' 'AMUA100-10-200_ERPAC_attign' '_' info.(currsubj).lohipitch{iset} '-pitch' '.mat'];
    save(fn,'attamuapac','ignamuapac','atttrials','igntrials')
    %%%%%
  end
  
  
  %%%%% get PAC for pitch
  % determine which stimuli were lo/hi
  switch info.(currsubj).leftrightelec
    case 'right'
      lotrials = find(attBRlostd);
      hitrials = find(attBRhistd);
      
    case 'left'
      lotrials = find(attBLlostd);
      hitrials = find(attBLhistd);
      
    otherwise
      error('woopsie')
  end
  
  % get fourfreq
  cfg = [];
  cfg.channel    = 'all';
  cfg.trials     = [lotrials' hitrials'];
  cfg.keeptrials = 'yes';
  cfg.pad        = 2; % pad out to 1 second;
  cfg.method     = 'mtmconvol';
  cfg.output     = 'fourier';
  cfg.toi        = toi; 
  cfg.foi        = 8:30;
  cfg.t_ftimwin  = 3./cfg.foi;
  cfg.taper      = 'hanning';
  freqphas = ft_freqanalysis(cfg,data);
  
  % get amua array
  amua = filtdata.trial([lotrials' hitrials']);
  amua = permute(cat(3,amua{:}),[3 1 2]);
  amua = cat(4,amua,amua);
  amua = permute(amua,[1 2 4 3]);
  
  % create freqdata
  freqdata = rmfield(freqphas,{'fourierspctrm','freq'});
  freqdata.phasfreqfour = freqphas.fourierspctrm;
  freqdata.phasfreq     = freqphas.freq;
  freqdata.ampfreqenv   = amua;
  freqdata.ampfreq      = [1 2]; % include both to be sure, in case of silly n=1 errors for some of the dimensions in roe_calc...
  freqdata.dimord       = 'rpt_chan_freq_time';
  clear freqphas
  
  % pre remove NaNs
  nanampenv    = isnan(squeeze(mean(mean(mean(freqdata.ampfreqenv,1),2),3)));
  nanfourierts = isnan(squeeze(mean(mean(mean(freqdata.phasfreqfour,1),2),3)));
  nonnantimepoints = logical(~nanampenv .* ~nanfourierts);
  freqdata.phasfreqfour = freqdata.phasfreqfour(:,:,:,nonnantimepoints);
  freqdata.ampfreqenv   = freqdata.ampfreqenv(:,:,:,nonnantimepoints);
  display(['new # time-points = ' num2str(sum(nonnantimepoints)) ' (removed ' num2str(numel(toi) - sum(nonnantimepoints)) ' timepoints)'])
  % save time axis to put back later
  timeinfo         = freqdata.time(nonnantimepoints);
  
  % split it
  % for lo
  freqlo = freqdata;
  freqlo.phasfreqfour = permute(freqlo.phasfreqfour(1:numel(lotrials),:,:,:), [4 2 3 1]);
  freqlo.ampfreqenv   = permute(freqlo.ampfreqenv(1:numel(lotrials),:,:,:), [4 2 3 1]);
  freqlo.time         = lotrials;
  % for ign
  freqhi = freqdata;
  freqhi.phasfreqfour = permute(freqhi.phasfreqfour(numel(lotrials)+1:end,:,:,:), [4 2 3 1]);
  freqhi.ampfreqenv   = permute(freqhi.ampfreqenv(numel(lotrials)+1:end,:,:,:), [4 2 3 1]);
  freqhi.time         = hitrials;
  clear freqdata
  
  % calc freqwphase and insert timeinfo
  cfg = [];
  cfg.handlenans    = 'conservative';
  cfg.normalization = 'withintrials';
  cfg.phasparam     = 'phasfreqfour';
  cfg.ampparam      = 'ampfreqenv';
  cfg.keeptrials    = 'yes';
  % for lo
  loamuapac = roe_calcwphaselocking(cfg,freqlo);
  loamuapac.timeinfo = timeinfo;
  % for hi
  hiamuapac = roe_calcwphaselocking(cfg,freqhi);
  hiamuapac.timeinfo = timeinfo;
  
  % save shit
  fn = [savepath currsubj '_' 'AMUA100-10-200_ERPAC_pitch' '.mat'];
  save(fn,'loamuapac','hiamuapac','lotrials','hitrials')
  %%%%%
  
  
end









% % DECOMP parafac
% for    isubj = 1    :nsubj
%   
%   % set
%   currsubj = info.subj{isubj};
%   disp(['working on ' currsubj])
%   
%   % go for it
%   currsubj = info.subj{isubj};
%   disp(['working on ' currsubj])
%   
%   % set global nwaycomp settings
%   cfg = [];
%   cfg.numitt             = 1000;
%   cfg.convcrit           = 1e-8;
%   cfg.randstart          = 10;
%   cfg.ncompestrandstart  = 10;
%   cfg.model              = 'parafac';
%   cfg.complexdims        = [1 1 0 0];
%   cfg.datparam           = 'wphaselockspctrm';
%   cfg.degencrit          = .7;
%   cfg.ncompest           = 'splithalf';
%   cfg.ncompestshcritval  = .7;
%   cfg.ncompeststart      = 1;
%   cfg.ncompestend        = 15;
%   cfg.ncompeststep       = 2;
%   cfg.datparam           = 'wplf';
%   cfg.ncompestshdatparam = 'wplfpart';
%   
%   % set fn part
%   nwayfnadd = ['_' 'nwaycomp' '_' cfg.ncompest '_' num2str(cfg.degencrit) '_' 'convcrit' num2str(cfg.ncompestshcritval)];
%   
%   % attention
%   for iset = 1:numel(info.(currsubj).lohipitch)
%     fn = [savepath currsubj '_' 'AMUA100-10-200_PAC_attign' '_' info.(currsubj).lohipitch{iset} '-pitch' '.mat'];
%     load(fn)
%     nwayfn = [savepath currsubj '_' 'AMUA100-10-200_PAC_attign' '_' info.(currsubj).lohipitch{iset} '-pitch' nwayfnadd '.mat'];
%     if ~exist(nwayfn,'file')
% 
%       % build freqdata
%       freqdata = [];
%       freqdata.wplf        = attamuapac.wphaselockspctrm + ignamuapac.wphaselockspctrm;
%       freqdata.wplfpart{1} = attamuapac.wphaselockspctrm;
%       freqdata.wplfpart{2} = ignamuapac.wphaselockspctrm;
%       freqdata.phasfreq    = attamuapac.phasfreq;
%       freqdata.ampfreq     = attamuapac.ampfreq;
%       freqdata.label       = attamuapac.label;
%       freqdata.dimord      = attamuapac.dimord;
%       freqdata.cfg         = [];
%  
%       % get it 
%       nwaycomp = nd_nwaydecomposition(cfg,freqdata);
%       % save
%       save(nwayfn,'nwaycomp')
%     end
%   end
%   
%   % pitch
%   fn = [savepath currsubj '_' 'AMUA100-10-200_PAC_pitch' '.mat'];
%   load(fn)
%   nwayfn = [savepath currsubj '_' 'AMUA100-10-200_PAC_pitch' nwayfnadd '.mat'];
%   if ~exist(nwayfn,'file') 
%     
%     % build freqdata
%     freqdata = [];
%     freqdata.wplf        = loamuapac.wphaselockspctrm + hiamuapac.wphaselockspctrm;
%     freqdata.wplfpart{1} = loamuapac.wphaselockspctrm;
%     freqdata.wplfpart{2} = hiamuapac.wphaselockspctrm;
%     freqdata.phasfreq    = loamuapac.phasfreq;
%     freqdata.ampfreq     = loamuapac.ampfreq;
%     freqdata.label       = loamuapac.label;
%     freqdata.dimord      = loamuapac.dimord;
%     freqdata.cfg         = [];
%     
%     % get it
%     nwaycomp = nd_nwaydecomposition(cfg,freqdata);
%     % save
%     save(nwayfn,'nwaycomp')
%   end
% end








function playground



%%%%%%%%%%%% BROWSE
% get details
info = rmr_abcaudecog_info;

% set stuff
savepath   = [info.savepath 'amua/'];

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % go for it
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  % attention
  for iset = 1%:numel(info.(currsubj).lohipitch)
    fn = [savepath currsubj '_' 'AMUA100-10-200_ERPAC_attign' '_' info.(currsubj).lohipitch{iset} '-pitch' '.mat'];
    load(fn)  
    
    % remove channels from lay
    load([info.savepath currsubj '_fieldtrip_layout.mat']);
    selchan = match_str(lay.label,attamuapac.label_old);
    lay.pos    = lay.pos(selchan,:);
    lay.label  = lay.label(selchan);
    lay.width  = lay.width(selchan);
    lay.height = lay.height(selchan);
    
    % hack for browsing
    attamuapac.wphaselockspctrm = permute(attamuapac.wphaselockspctrm(:,:,:,1,:),[2 3 5 1 4]); 
    attamuapac.dimord   = 'ampchan_phaschan_ampfreq_phasfreq';
    attamuapac.ampfreq  = attamuapac.phasfreq;
    attamuapac.phasfreq = attamuapac.timeinfo;
    ignamuapac.wphaselockspctrm = permute(ignamuapac.wphaselockspctrm(:,:,:,1,:),[2 3 5 1 4]);
    ignamuapac.dimord = 'ampchan_phaschan_ampfreq_phasfreq';
    ignamuapac.ampfreq  = ignamuapac.phasfreq;
    ignamuapac.phasfreq = ignamuapac.timeinfo;
    
    % att
    cfg = [];
    cfg.layout     = lay;
    cfg.windowname = [currsubj '  att' '_' info.(currsubj).lohipitch{iset}];
    roe_browsefreqwphase(cfg,attamuapac)
    % ign
    cfg = [];
    cfg.layout     = lay;
    cfg.windowname = [currsubj '  ign' '_' info.(currsubj).lohipitch{iset}];
    roe_browsefreqwphase(cfg,ignamuapac)
  end
  
  % [itch
  fn = [savepath currsubj '_' 'AMUA100-10-200_ERPAC_pitch' '.mat'];
  load(fn)
  
  % remove channels from lay
  load([info.savepath currsubj '_fieldtrip_layout.mat']);
  selchan = match_str(lay.label,freqwphase.label_old);
  lay.pos    = lay.pos(selchan,:);
  lay.label  = lay.label(selchan);
  lay.width  = lay.width(selchan);
  lay.height = lay.height(selchan);
  
  % calculate shuffle percentile 0.95 and set everything below to 0/NaN
  % for lo
  shufmean = squeeze(mean(abs(loamuapac.wphaselockspctrm(1:end-1,:,:,:,:)),1));
  shufstd  = squeeze(std(abs(loamuapac.wphaselockspctrm(1:end-1,:,:,:,:)),[],1));
  p95 = norminv(0.95,shufmean,shufstd);
  loamuapac.wphaselockspctrm = squeeze(loamuapac.wphaselockspctrm(end,:,:,:,:));
  loamuapac.wphaselockspctrm(loamuapac.wphaselockspctrm<p95) = 0;
  % for ign
  shufmean = squeeze(mean(abs(hiamuapac.wphaselockspctrm(1:end-1,:,:,:,:)),1));
  shufstd  = squeeze(std(abs(hiamuapac.wphaselockspctrm(1:end-1,:,:,:,:)),[],1));
  p95 = norminv(0.95,shufmean,shufstd);
  hiamuapac.wphaselockspctrm = squeeze(hiamuapac.wphaselockspctrm(end,:,:,:,:));
  hiamuapac.wphaselockspctrm(hiamuapac.wphaselockspctrm<p95) = 0;
   
  % lo
  cfg = [];
  cfg.layout     = lay;
  cfg.windowname = [currsubj '  lo'];
  roe_browsefreqwphase(cfg,loamuapac)
  % hi
  cfg = [];
  cfg.layout     = lay;
  cfg.windowname = [currsubj '  hi'];
  roe_browsefreqwphase(cfg,hiamuapac)
end

























%%%%%%%%%%%% PLOT WITHIN CHAN
% get details
info = rmr_abcaudecog_info;

% set stuff
savepath   = [info.savepath 'amua/'];

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % go for it
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  load([info.savepath currsubj '_fieldtrip_layout.mat']);

  % attention
  for iset = 1:numel(info.(currsubj).lohipitch)
    fn = [savepath currsubj '_' 'AMUA100-10-200_PAC_attign' '_' info.(currsubj).lohipitch{1} '-pitch' '.mat'];
    load(fn)
       
    % att
    figure('name',[currsubj ' att' '_' info.(currsubj).lohipitch{iset}],'numbertitle','off')
    hold on
    % draw outlines
    for ioutline = 1:numel(lay.outline)
      xline = lay.outline{ioutline}(:,1);
      yline = lay.outline{ioutline}(:,2);
      line(xline,yline,'linewidth',1,'linestyle','--')
    end
    % draw data
    label = freqwphase.label_old;
    nchan = numel(label);
    for ichan = 1:nchan
      % get local axis
      layind = strcmp(lay.label,label{ichan});
      hpos   = lay.pos(layind,1);
      vpos   = lay.pos(layind,2);
      width  = lay.width(layind);
      height = lay.height(layind);
      
      % plot patch
      C = abs(squeeze(attamuapac.wphaselockspctrm(ichan,ichan,:,:)));
      clim = [0 0.05];
      ft_plot_matrix(C,'hpos',hpos,'vpos',vpos,'height',height,'width',width,'clim',clim)
    end
    axis off
    
    % ign
    figure('name',[currsubj ' ign' '_' info.(currsubj).lohipitch{iset}],'numbertitle','off')
    hold on
    % draw outlines
    for ioutline = 1:numel(lay.outline)
      xline = lay.outline{ioutline}(:,1);
      yline = lay.outline{ioutline}(:,2);
      line(xline,yline,'linewidth',1,'linestyle','--')
    end
    % draw data
    label = freqwphase.label_old;
    nchan = numel(label);
    for ichan = 1:nchan
      % get local axis
      layind = strcmp(lay.label,label{ichan});
      hpos   = lay.pos(layind,1);
      vpos   = lay.pos(layind,2);
      width  = lay.width(layind);
      height = lay.height(layind);
      
      % plot patch
      C = abs(squeeze(ignamuapac.wphaselockspctrm(ichan,ichan,:,:)));
      clim = [0 0.05];
      ft_plot_matrix(C,'hpos',hpos,'vpos',vpos,'height',height,'width',width,'clim',clim)
    end
    axis off
    
  end
  
  % [itch
  fn = [savepath currsubj '_' 'AMUA100-10-200_PAC_pitch' '.mat'];
  load(fn)
  
    % lo
    figure('name',[currsubj ' lo' '_' info.(currsubj).lohipitch{iset}],'numbertitle','off')
    hold on
    % draw outlines
    for ioutline = 1:numel(lay.outline)
      xline = lay.outline{ioutline}(:,1);
      yline = lay.outline{ioutline}(:,2);
      line(xline,yline,'linewidth',1,'linestyle','--')
    end
    % draw data
    label = freqwphase.label_old;
    nchan = numel(label);
    for ichan = 1:nchan
      % get local axis
      layind = strcmp(lay.label,label{ichan});
      hpos   = lay.pos(layind,1);
      vpos   = lay.pos(layind,2);
      width  = lay.width(layind);
      height = lay.height(layind);
      
      % plot patch
      C = abs(squeeze(loamuapac.wphaselockspctrm(ichan,ichan,:,:)));
      clim = [0 0.05];
      ft_plot_matrix(C,'hpos',hpos,'vpos',vpos,'height',height,'width',width,'clim',clim)
    end
    axis off
    
    % hi
    figure('name',[currsubj ' hi' '_' info.(currsubj).lohipitch{iset}],'numbertitle','off')
    hold on
    % draw outlines
    for ioutline = 1:numel(lay.outline)
      xline = lay.outline{ioutline}(:,1);
      yline = lay.outline{ioutline}(:,2);
      line(xline,yline,'linewidth',1,'linestyle','--')
    end
    % draw data
    label = freqwphase.label_old;
    nchan = numel(label);
    for ichan = 1:nchan
      % get local axis
      layind = strcmp(lay.label,label{ichan});
      hpos   = lay.pos(layind,1);
      vpos   = lay.pos(layind,2);
      width  = lay.width(layind);
      height = lay.height(layind);
      
      % plot patch
      C = abs(squeeze(hiamuapac.wphaselockspctrm(ichan,ichan,:,:)));
      clim = [0 0.05];
      ft_plot_matrix(C,'hpos',hpos,'vpos',vpos,'height',height,'width',width,'clim',clim)
    end
    axis off
    
end
















