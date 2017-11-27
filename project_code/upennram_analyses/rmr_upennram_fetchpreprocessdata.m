function data = rmr_upennram_fetchpreprocessdata(cfgprepro,currsubj,info)

% set defaults
cfgprepro.detrend        = ft_getopt(cfgprepro, 'detrend',      'no');
cfgprepro.demean         = ft_getopt(cfgprepro, 'demean',       'yes');
cfgprepro.useenc         = ft_getopt(cfgprepro, 'useenc',       true);
cfgprepro.userec         = ft_getopt(cfgprepro, 'userec',       true);
cfgprepro.encduration    = ft_getopt(cfgprepro, 'encduration',  1.6);
cfgprepro.encprestim     = ft_getopt(cfgprepro, 'encprestim',   0);
cfgprepro.encpoststim    = ft_getopt(cfgprepro, 'encpoststim',  0);
cfgprepro.recduration    = ft_getopt(cfgprepro, 'recduration',  0.5);
cfgprepro.recprestim     = ft_getopt(cfgprepro, 'recprestim',   0);
cfgprepro.recpoststim    = ft_getopt(cfgprepro, 'recpoststim',  0);
cfgprepro.reref          = ft_getopt(cfgprepro, 'reref',        'yes');
cfgprepro.bsfilter       = ft_getopt(cfgprepro, 'bsfilter',     'yes');
cfgprepro.lpfilter       = ft_getopt(cfgprepro, 'lpfilter',     'yes');


% combine data over sessions and preprocess
datacmb = [];
imont   = 1; % HARD-CODED, currently no other montages present (and should be treated largely as a different subject anyway)
for   isess = 1  :numel(info.sessions.(currsubj){imont})
  disp(['working on ' currsubj ' session ' num2str(isess) ' of ' num2str(numel(info.sessions.(currsubj){imont}))])
  
  % fetch header
  hdr   = read_upennram_header([info.datapath info.sessions.(currsubj){imont}(isess).headerfile]);
  
  % obtain segmentation
  cfg = []; % start with an empty cfg
  cfg.header      = hdr;
  cfg.event       = read_upennram_event([info.datapath info.sessions.(currsubj){imont}(isess).eventfile]);
  cfg.encduration = cfgprepro.encduration; % during encoding, the period, in seconds, after/before pre/poststim periods
  cfg.recduration = cfgprepro.recduration; % during   recall, the period, in seconds, after/before pre/poststim periods
  cfg.encprestim  = cfgprepro.encprestim;  % during encoding, the period, in seconds, before word onset that is additionally cut out (t=0 will remain word onset)
  cfg.encpoststim = cfgprepro.encpoststim; % during encoding, the period, in seconds, after cfg.encduration, that is additionally cut out
  cfg.recprestim  = cfgprepro.recprestim;  % during   recall, the period, in seconds, before word onset that is additionally cut out (t=0 will remain word onset)
  cfg.recpoststim = cfgprepro.recpoststim; % during   recall, the period, in seconds, after cfg.recduration, that is additionally cut out
  trl = rmr_upennram_trialfun(cfg); % obtain the trl matrix, which contains the segmentation details (Note, this function can also be called from with ft_definetrial)
  % keep only parts requested
  if ~cfgprepro.userec
    remind = trl(:,4)==2;
    trl(remind,:) = [];
  end
  if ~cfgprepro.useenc
    remind = trl(:,4)==1;
    trl(remind,:) = [];
  end
  
  % reject artifacts using a trick and ft_rejectartifact
  cfg = [];
  cfg.trl          = trl;
  cfg.headerfile   = [info.datapath info.sessions.(currsubj){imont}(isess).headerfile];
  cfg.headerformat = 'read_upennram_header';
  cfg.artfctdef.visual.artifact =  info.(currsubj).artfctdef.visual.(hdr.channelfile{1}(1:end-4));
  cfgout = ft_rejectartifact(cfg);
  trl = cfgout.trl; 

  % read in data and preprocess
  cfg = [];
  cfg.trl           = trl;
  cfg.datafile      = [info.datapath info.sessions.(currsubj){imont}(isess).datadir];
  cfg.dataformat    = 'read_upennram_data';
  cfg.headerfile    = [info.datapath info.sessions.(currsubj){imont}(isess).headerfile];
  cfg.headerformat  = 'read_upennram_header';
  cfg.padding       = (max(cfg.trl(:,2)-cfg.trl(:,1))+1)./hdr.Fs + max((2*info.(currsubj).bsfilt.edgeartlen),2*0.2);
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
  cfg.lpfreq        = 300;
  cfg.lpfilttype    = 'firws';
  cfg.lpfiltdir     = 'onepass-zerophase';
  cfg.lpfiltwintype = 'blackman';
  cfg.lpfiltord     = round(0.2*hdr.Fs);
  cfg.usefftfilt    = 'yes';
  % set the channels to be used
  cfg.channel       = setdiff(hdr.label,ft_channelselection(info.(currsubj).badchan,hdr.label));
  % reref
  if strcmp(cfgprepro.reref,'yes')
    cfg.reref       = 'yes';
    cfg.refchannel  = cfg.channel;
  end
  % analysis specific options
  cfg.demean          = cfgprepro.demean;
  cfg.detrend         = cfgprepro.detrend;
  %%%
  % allow for overide of some of the default settings by cfgprepro
  fieldstocopy = setdiff(fieldnames(cfgprepro),{'useenc','userec','encduration','encprestim','encpoststim','recduration','recprestim','recpoststim','detrend','demean','reref','bsfilter','lpfilter'});
  cfg = copyfields(cfgprepro,cfg,fieldstocopy);
  %%%
  % make it so
  data = ft_preprocessing(cfg);
  
  % keep data for combining
  datacmb{isess} = data;
  clear data
end % isess

% gather fsample
fsample = cell(size(datacmb));
for isess = 1:numel(datacmb)
  fsample{isess} = datacmb{isess}.fsample;
end

% combine sessions
if numel(datacmb)>1
  data = ft_appenddata([],datacmb{:});
  clear datacmb
else
  data = datacmb{1};
  clear datacmb
end

% put back fsample
if numel(fsample)>1 && isequal(fsample{:})
  data.fsample = fsample{1};
end









