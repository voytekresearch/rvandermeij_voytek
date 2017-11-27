function data = rmr_predfaceval_fetchpreprocessdata(cfgprepro,currsubj,info)

% set defaults
cfgprepro.detrend     = ft_getopt(cfgprepro, 'detrend',      'no');
cfgprepro.demean      = ft_getopt(cfgprepro, 'demean',       'yes');
cfgprepro.prestim     = ft_getopt(cfgprepro, 'prestim',      0);
cfgprepro.poststim    = ft_getopt(cfgprepro, 'poststim',     0);
cfgprepro.reref       = ft_getopt(cfgprepro, 'reref',        'no');
cfgprepro.refchannel  = ft_getopt(cfgprepro, 'refchannel',   []);
cfgprepro.bsfilter    = ft_getopt(cfgprepro, 'bsfilter',     'yes');
cfgprepro.lpfilter    = ft_getopt(cfgprepro, 'lpfilter',     'yes');


% combine data over sessions and preprocess
datacmb = [];
for   isess = 1  :numel(info.sessiondata.(currsubj))
  currevents = info.sessionevents.(currsubj){isess};
  currdata   = info.sessiondata.(currsubj){isess};
  disp(['working on ' currsubj ': ' 'data = ' currdata ', events = ' currevents])
  
  % obtain segmentation
  cfg = [];
  cfg.datafile  = [info.datapath currsubj '/' currdata];
  cfg.eventfile = [info.datapath currsubj '/' currevents];
  cfg.diodechan = info.sessiondiode.(currsubj){isess};
  cfg.threshold = info.sessionthresh.(currsubj){isess};
  cfg.prestim   = cfgprepro.prestim;     % the period, in seconds, before CUE ONSET that is additionally cut out (t=0 will remain CUE ONSET)
  cfg.poststim  = cfgprepro.poststim;    % the period, in seconds, after FACE ONSET that is additionally cut out (FACE ONSET is kept in data.trialinfo, see above)
  cfg.debugflg  = false;
  trl = rmr_predfaceval_definetrials(cfg); % obtain the trl matrix, which contains the segmentation details
  
  % reject using a trick and ft_rejectartifact
  cfg = [];
  cfg.trl          = trl;
  cfg.headerfile   = [info.datapath currsubj '/' currdata];
  [dum, main, ext] = fileparts(currdata);
  cfg.artfctdef.roemer.artifact = info.(currsubj).artfctdef.visual.(([ext(2:end) main]));
  cfgout = ft_rejectartifact(cfg);
  trl = cfgout.trl;
  
  % remove trials with cuefaceinterval above 1 second (next step is 1.5)
  remind = (trl(:,8) - trl(:,11)) > 1.25;
  trl(remind,:) = [];
  % remove trials with incorrect response
  remind = trl(:,6) == 0;
  trl(remind,:) = [];
  
  % readin data
  % get header for channel detection
  hdr = ft_read_header([info.datapath currsubj '/' currdata]);
  % read in and filter
  cfg = [];
  cfg.trl           = trl;
  cfg.datafile      = [info.datapath currsubj '/' currdata];
  cfg.padding       = (max(cfg.trl(:,2)-cfg.trl(:,1))+1)./hdr.Fs + max((2*info.(currsubj).bsfilt.edgeartlen),2*0.2);
  cfg.reref         = cfgprepro.reref;
  cfg.refchannel    = cfgprepro.refchannel;
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
  cfg.lpfreq        = 540;
  cfg.lpfilttype    = 'firws';
  cfg.lpfiltdir     = 'onepass-zerophase';
  cfg.lpfiltwintype = 'blackman';
  cfg.lpfiltord     = round(0.2*hdr.Fs);
  cfg.usefftfilt    = 'yes';
  % set the channels to be removed
  if isempty(strfind(currdata,'edf'))
    cfg.channel       = setdiff(hdr.label,info.(currsubj).badchan);
  else
    cfg.channel       = setdiff(hdr.label,rmr_predfaceval_besatoedfchannelnames(info.(currsubj).badchan));
  end
  % analysis specific options
  cfg.demean          = cfgprepro.demean;
  cfg.detrend         = cfgprepro.detrend;
  %%%
  % allow for overide of some of the default settings by cfgprepro
  fieldstocopy = setdiff(fieldnames(cfgprepro),{'prestim','poststim'});
  cfg = copyfields(cfgprepro,cfg,fieldstocopy);
  %%%
  % make it so
  data = ft_preprocessing(cfg);
  % clean channel names if necessary
  data.label = rmr_predfaceval_edftobesachannelnames(data.label);
  
  % keep data for combining
  datacmb{isess} = data;
end % isess

% combine sessions
if numel(datacmb)>1
  data = ft_appenddata([],datacmb{:});
  clear datacmb
else
  data = datacmb{1};
  clear datacmb
end

% remove trials with above x*mad +/- median
rt     = data.trialinfo(:,4);
datmed = median(rt);
datmad = mad(rt,1);
remind = (rt < (datmed-6*datmad)) | (rt > (datmed+6*datmad));
cfg = [];
cfg.trial = find(~remind);
data = ft_preprocessing(cfg,data);








