function data = rmr_abcaudecog_fetchpreprocessdata(cfgprepro,currsubj,info)

% set defaults
cfgprepro.detrend     = ft_getopt(cfgprepro, 'detrend',      'no');
cfgprepro.demean      = ft_getopt(cfgprepro, 'demean',       'yes');
cfgprepro.prestim     = ft_getopt(cfgprepro, 'prestim',      0);
cfgprepro.poststim    = ft_getopt(cfgprepro, 'poststim',     0);
cfgprepro.reref       = ft_getopt(cfgprepro, 'reref',        'no');
cfgprepro.refchannel  = ft_getopt(cfgprepro, 'refchannel',   []);
cfgprepro.bsfilter    = ft_getopt(cfgprepro, 'bsfilter',     'yes');
cfgprepro.lpfilter    = ft_getopt(cfgprepro, 'lpfilter',     'yes');


% comb sessions and preprocess
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
  % set the bandstops
  cfg.bsfilter      = cfgprepro.bsfilter;
  bsfiltpeak  = info.(currsubj).bsfilt.peak;
  bsfiltbandw = info.(currsubj).bsfilt.halfbandw;
  cfg.bsfreq        = [bsfiltpeak-bsfiltbandw;  bsfiltpeak+bsfiltbandw]';
  cfg.bsfilttype    = 'but';
  cfg.bsfiltdir     = 'twopass';
  cfg.bsfiltord     = 2;
  % set the low pass
  cfg.lpfilter      = cfgprepro.lpfilter;
  if data.fsample==1000
    cfg.lpfreq     = 350;
  else
    cfg.lpfreq     = 540;
  end
  cfg.lpfilttype    = 'firws';
  cfg.lpfiltdir     = 'onepass-zerophase';
  cfg.lpfiltwintype = 'blackman';
  cfg.lpfiltord     = round(0.2*data.fsample);
  cfg.usefftfilt    = 'yes';
  % analysis specific options
  cfg.demean     = cfgprepro.demean;
  cfg.detrend    = cfgprepro.detrend;
  cfg.reref      = cfgprepro.reref;
  cfg.refchannel = cfgprepro.refchannel;
  cfg.channel    = setdiff(data.label,info.(currsubj).badchan);
  %%%
  % allow for overide of some of the default settings by cfgprepro
  fieldstocopy = setdiff(fieldnames(cfgprepro),{'prestim','poststim'});
  cfg = copyfields(cfgprepro,cfg,fieldstocopy);
  %%%
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
  %%%%%%%%%%%
  % add presetim
  trl(:,1) = trl(:,1) - round(cfgprepro.prestim .* data.fsample);
  trl(:,3) = -round(cfgprepro.prestim .* data.fsample);
  % post stim
  trl(:,2) = trl(:,2) + round(cfgprepro.poststim .* data.fsample);
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
  cfg = [];
  cfg.artfctdef.visual.artifact = info.(currsubj).artfctdef.visual.(currsess);
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