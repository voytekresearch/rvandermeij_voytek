function wplfdata = rmr_phaseamplitudecoupling(cfg,data)

% RMR_PHASEAMPLITUDECOUPLING calculates amplitude weigthed Phase-Locking Factors
% that quantify coupling between the amplitude envelope of one frequency and the 
% phase (and amplitude) of another. In the following, amplitude(-providing) frequencies 
% refer to the former, and phase(-providing) frequencies refer to the latter.
% Spectral analysis is performed inside this function using ft_specest_mtmconvol 
% (from the FieldTrip SPECEST module).
%
% Use as
%   [freqwphase] = rmr_phaseamplitudecoupling(cfg, data)
%
% The input data should be organised in a structure as obtained from FT_PREPROCESSING. 
% The configuration structure contains options as specified below.
%
%        Input options:
%          - Spectral Analysis options - 
%    cfg.ampfreq           = 1xNampfreq vector, amplitude frequencies to be estimated in Hz
%    cfg.phasfreq          = 1xNphasfreq vector, phase frequencies to be estimated in Hz
%    cfg.ampfreqtimwin     = 1xNampfreq vector, time-window to use for estimating amplitude frequencies (analogous to t_ftimwin)
%    cfg.phasfreqtimwin    = 1xNphasfreq vector, time-window to use for estimating amplitude frequencies (analogous to t_ftimwin)
%    cfg.ampfreqtaper      = string, function to taper time-window to be used in spectral analysis with (default = hanning)
%    cfg.phasfreqtaper     = string, function to taper time-window to be used in spectral analysis with (default = hanning)
%    cfg.amptstohfats      = 'yes', 'no' (default = 'no'), use above settings to combine amplitude frequencies into one HFA time series
%                                  HFA will be estimated by averaging (over frequencies) z-scored amplitude time-series. Afterwards, wplfdata.ampfreq will be set to the average frequency.
%       Note1: data will be zero-padded to an integer number of seconds during spectral analysis 
%       Note2: the maximal amount of data of each trial will be used, data selection needs to occur prior to calling rmr_phaseamplitudecoupling 
%    
%          - General options - 
%    cfg.normalization     = 'withintrials', 'overtrials', 'none' 
%                                  withintrials: wPLFs are normalized per trial, and wPLFs in output are the average over trials (if keeptrials = no)
%                                  overtrials:  wPLFs are normalized using the SSQ of all trials combined, i.e. as if all trials were a single one
%                                  none: wPLFs are the raw cross-products of the amplitude envelopes and phases (and amplitudes), summed over trials
%    cfg.handlespecedges   = 'conservative' or 'liberal', handling of time points at edges of trials where wavelets are not fully immersed in the data
%                                  convervative: for every frequency-pair, only time-points are taken for which the wavelet with longest time-domain length was fully immersed in the data
%                                  liberal: for every frequency-pair, the maximal amount of time points are used
%    cfg.keeptrials        = 'yes' or 'no' (default = 'no')
%                                   Note: if normalization = 'overtrials' and keeptrials = 'yes', than trial-specific wPLFs will be normalized as if
%                                   keeptrials = 'no', and the "final" wPLFs that are bound between 0 and 1 are obtained by summing wPLFs over trials
%    cfg.centerampfreq     = 'no','yes' (default = 'yes'), center the amplitude frequency time-series prior to calculating wPLF
%                                   No:  wPLFs, between conditions, are sensitive to additive differences of the amplitude envelope, and insenstive to multiplicative differences
%                                   Yes: wPLFs, between conditions, are insensitive to both additive and multiplicative differences of the amplitude envelope 
%    cfg.phasfreqnoamp     = 'no','yes' (default = 'no'), only use the phase of phase frequencies and remove amplitude (default = 'no') 
%                                   Note: phasfreqnoamp = 'yes' combined with centerampfreq = 'no' will make wPLFs equivalent to "Ozkurt's" coefficient (see Ozkurt & Schnittzler 2011 JNeuroMethods)
%    cfg.nsurrbytrials     = empty, or number of surrogate wPLFs to obtain by calculating wPLFs between amplitude and phase frequencies of randomly different trials
%    cfg.nsurrbytime       = empty, or number of surrogate wPLFs to obtain by calculating wPLFs between randomly circular time-shifted amplitude and phase frequencies (max shift = length of time series)
%                                   Note: only of of the surrogate options can be used
%    cfg.trials            = 'all' or a selection given as a 1xN vector (default = 'all')
%    cfg.channel           = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details (all combinations will be calculated)
%    cfg.precision         = 'single' or 'double', precision of wPLFs (computations will be done at the precision of the raw data)
%    cfg.randomseed        = scalar, seed of the random number generator
%
%
%        Output:
%                 label: 1xNchan celly-array, channels labels of used channels
%               ampfreq: 1xNampfreq vector, amplitude frequencies in Hz
%              phasfreq: 1xNphasfreq vector, phase frequencies in Hz
%                dimord: identity of dimensions of wplf
%                  wplf: complex-valued array of wPLFs
%         surrogatewplf: 1xNsurrogates cell-array, each containing an array of surrogate wPLFs
%                   dof: NtrialxNampfreqxNphasfreq array, number time-points for each of the frequency pairs used of each trial
%             trialinfo: trialinfo field from input data for trials that were used
%                   cfg: cfg used and its history
%
%
%
%               To facilitate data-handling and distributed computing you can use
%    cfg.inputfile         =  ...
%    cfg.outputfile        =  ...
%               If you specify one of these (or both) the input data will be read from a *.mat
%               file on disk and/or the output data will be written to a *.mat file. These mat
%               files should contain only a single variable, corresponding with the
%               input/output structure.
%
%
% This function depends on the FieldTrip toolbox, which can be found at www.fieldtriptoolbox.org.
%
% Copyright (C) 2016-present, Roemer van der Meij
%


% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug 
ft_preamble loadvar data 

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check for required options
ft_checkconfig(cfg,'required',{'ampfreq','phasfreq','ampfreqtimwin','phasfreqtimwin'});
ft_checkconfig(cfg,'required',{'normalization','handlespecedges'});

% set default options
cfg.normalization     = ft_getopt(cfg, 'normalization',   []);
cfg.handlespecedges   = ft_getopt(cfg, 'handlespecedges',      []);
cfg.ampfreqtaper      = ft_getopt(cfg, 'ampfreqtaper',    'hanning');
cfg.phasfreqtaper     = ft_getopt(cfg, 'phasfreqtaper',   'hanning');
cfg.amptstohfats      = ft_getopt(cfg, 'amptstohfats',    'no');
cfg.keeptrials        = ft_getopt(cfg, 'keeptrials',      'no');
cfg.nsurrbytrials     = ft_getopt(cfg, 'nsurrbytrials',     []);
cfg.nsurrbytime       = ft_getopt(cfg, 'nsurrbytime',       []);
cfg.centerampfreq     = ft_getopt(cfg, 'centerampfreq',    'yes');
cfg.phasfreqnoamp     = ft_getopt(cfg, 'phasfreqnoamp',    'no');
cfg.trials            = ft_getopt(cfg, 'trials',          'all');
cfg.channel           = ft_getopt(cfg, 'channel',         'all');
cfg.precision         = ft_getopt(cfg, 'precision',       'double');
cfg.randomseed        = ft_getopt(cfg, 'randomseed',      sum(clock.*1e6));

% check (depedendent) options
if (numel(cfg.ampfreq) ~= numel(cfg.ampfreqtimwin)) || (numel(cfg.phasfreq) ~= numel(cfg.phasfreqtimwin))
  error('cfg.amp/phasfreq and cfg.amp/phasfreqtimwin need to be the same size')
end
if ~isempty(cfg.nsurrbytrials) && ~isempty(cfg.nsurrbytime)
  error('only one of the surrogate options can be used at the same time')
end
if (~isempty(cfg.nsurrbytrials) && ~isnumeric(cfg.nsurrbytrials)) || (~isempty(cfg.nsurrbytime) && ~isnumeric(cfg.nsurrbytime))
  error('surrogate specification needs to be a single number')
end
if (~isempty(cfg.nsurrbytrials) || ~isempty(cfg.nsurrbytime)) && istrue(cfg.keeptrials)
  error('keeping trials while computing surrogate wPLFs is not supported (nor advised)')
end
if strcmp(cfg.ampfreqtaper,'dpss') || strcmp(cfg.phasfreqtaper,'dpss')
  error('multitapering is not supported')
end
if istrue(cfg.amptstohfats) && any(diff(cfg.ampfreqtimwin))
  error('the time-window for every amplitude frequency should be equal when calculting HFA')
end
if strcmp(cfg.normalization,'overtrials') && istrue(cfg.centerampfreq)
  error('mean centering amplitude frequencies is not supported when normalizing wPLFs over trials')
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw', 'raw+comp', 'mvar'}, 'hassampleinfo', 'yes');

% select trials and channels of interest
data = ft_selectdata(rmfield(cfg,setdiff(fieldnames(cfg),{'trials','channel'})), data);

% check options that cannot be used when trials do not have identical time axes
if strcmp(cfg.normalization,'overtrials') && ~isequal(data.time{:})
  error('normalizing wPLFs over trials is not supported with trials with different time axes')
end
if strcmp(cfg.handlespecedges,'conservative') && ~isequal(data.time{:})
  error('conservative handling of edges is not possible for trials with different time axes')
end
if ~isempty(cfg.nsurrbytrials) && ~isequal(data.time{:})
  error('computing surrogate wPLFs by shuffling trials is not possible when trials have different time axes')
end


% apply random seed
rng(cfg.randomseed);

% parse surrogate options
dosurr = (~isempty(cfg.nsurrbytrials) && isnumeric(cfg.nsurrbytrials) && cfg.nsurrbytrials>0) || (~isempty(cfg.nsurrbytime) && isnumeric(cfg.nsurrbytime) && cfg.nsurrbytime>0);


% provide some display (and also use as option check)
switch cfg.normalization  
  case 'withintrials'
    disp('normalizing wPLFs within trials')
  case 'overtrials'
    disp('normalizing wPLFs over trials')
  case 'none'
    disp('not normalizing wPLFs')
  otherwise
    error(['cfg.normalization = ' cfg.normalization ' is not supported'])
end
switch cfg.handlespecedges
  case 'liberal'
    disp('calculating wPLFs using maximally available number of time points')
  case 'conservative'
    disp('calculating wPLFs using only time points available to all frequency pairs')
  otherwise
    error(['cfg.handlespecedges = ' cfg.handlespecedges ' is not supported'])
end 


%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
% get wPLFs
% first get surrogate wPLFs
if dosurr
  nsurr = max([cfg.nsurrbytrials cfg.nsurrbytime]);
  
  % we have to do preallocation here, in order to prevent out of memory errors halfway through the shuffling
  wplfsize = [numel(data.label)  numel(data.label) numel(cfg.ampfreq) numel(cfg.phasfreq)];
  if istrue(cfg.keeptrials)
    wplfsize = [wplfsize numel(data.trial)];
  end
  surrogatewplf = cell(1,nsurr);
  for isurr = 1:nsurr
    surrogatewplf{isurr} = complex(NaN(wplfsize,cfg.precision),NaN(wplfsize,cfg.precision));
  end
  
  % get the surrogate wPLFs
  for isurr = 1:nsurr
    disp(['obtaining surrogate wPLFs set ' num2str(isurr) ' of ' num2str(nsurr)]);
    surrogatewplf{isurr} = calcwplfs(cfg,data);
  end
end
% get real wPLFs
disp('obtaining wPLFs');
[wplf,dof] = calcwplfs(cfg,data);
%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%




% create dimord
if istrue(cfg.keeptrials)
  dimord = 'rpt_ampchan_phaschan_ampfreq_phasfreq';
else
  dimord = 'ampchan_phaschan_ampfreq_phasfreq';
end

% create output
wplfdata = [];
wplfdata.label           = data.label;
if istrue(cfg.amptstohfats)
  wplfdata.ampfreq       = mean(cfg.ampfreq);
else
  wplfdata.ampfreq       = cfg.ampfreq;
end
wplfdata.phasfreq        = cfg.phasfreq;
wplfdata.dimord          = dimord;
wplfdata.wplf            = wplf;
if dosurr
  wplfdata.surrogatewplf = surrogatewplf;
end
wplfdata.dof             = dof;
wplfdata.trialinfo       = data.trialinfo;
  

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous data
ft_postamble history wplfdata
ft_postamble savevar wplfdata 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wplf,dof] = calcwplfs(cfg,data)
%
%
% Obtain wPLFs
%
%

% set sizes
nampfreq  = numel(cfg.ampfreq);
nphasfreq = numel(cfg.phasfreq);
nchan     = numel(data.label);
ntrial    = numel(data.trial);
if istrue(cfg.amptstohfats)
  nampfreq = 1;
end

% flg surrs
dosurrbytrials = ~isempty(cfg.nsurrbytrials); % accuracy of both is checked at the top of the main
dosurrbytime   = ~isempty(cfg.nsurrbytime);

% pre-allocate output variables
if ~istrue(cfg.keeptrials)
  wplfsize = [nchan nchan nampfreq nphasfreq];
else
  wplfsize = [ntrial nchan nchan nampfreq nphasfreq]; % ideally, speedwise, trial should be third dimension, but fuck it
end
wplf = complex(zeros(wplfsize,cfg.precision),zeros(wplfsize,cfg.precision));
dof  = zeros([ntrial nampfreq nphasfreq]);

% pre-allocate other variables
switch cfg.normalization
  case 'overtrials'
    ssqampts  = zeros(nchan,nampfreq,nphasfreq);
    ssqphasts = zeros(nchan,nampfreq,nphasfreq);
  case 'withintrials'
    % do nothing
end

% determine padding for spectral analysis
padding = ceil(max(cellfun(@numel,data.time)) ./ data.fsample);

% determine timeoi for spectral analysis
if strcmp(cfg.handlespecedges,'conservative')
  % all trials have been confirmed to have indentical timing above, so fine to take first
  timeoi  = data.time{1};
  timeboi = 1:numel(timeoi);
  % determine longest wavelet
  maxtimwin  = max(max(cfg.ampfreqtimwin),max(cfg.phasfreqtimwin));
  nsmptimiwn = maxtimwin .* data.fsample;
  % find time points that are fully immersed
  ind        = find((timeboi >=  (nsmptimiwn ./ 2)) & (timeboi < numel(timeoi) - (nsmptimiwn ./2)));
  timeoi     = timeoi(ind);
end

% obtain HFA if requested
if istrue(cfg.amptstohfats)
  % get hfa
  freqoi = cfg.ampfreq;
  timwin = cfg.ampfreqtimwin(1); % equality of timwin is checked at the top of main
  taper  = cfg.ampfreqtaper;
  if strcmp(cfg.handlespecedges,'conservative')
    hfatimeoi = timeoi;
  else
    hfatimeoi = 'all';
  end
  hfa = rmr_highfrequencyactivity(data,freqoi,hfatimeoi,timwin,taper);
end
  
  
% set trial indices to loop over (will be replaced when shuffling
trialind = [];
trialind(:,1) = 1:ntrial;
trialind(:,2) = 1:ntrial;

%%% Surrogate by trial shuffling %%%
if dosurrbytrials
  trialind(:,1) = randperm(ntrial);
  trialind(:,2) = randperm(ntrial);
end
%%% Surrogate by trial shuffling %%%

% calculate wPLFs
for itrial = 1:ntrial
  disp(['calculating wPLFs for trial ' num2str(itrial) ' of ' num2str(ntrial)]);
 
  % set trials to use
  trialsel1 = trialind(itrial,1);
  trialsel2 = trialind(itrial,1);
  
  % obtain Fourier coefficients for amplitude and phase frequencies
  opt = {'pad',padding,'dimord','chan_time_freqtap','verbose',0};
  if strcmp(cfg.handlespecedges,'conservative')
    opt = [opt {'timeoi',timeoi}];
  else
    opt = [opt {'timeoi',data.time{trialsel1}}]; % note, trial shuffling is only allowed when all trials have identical time axes, so it is fine that the same timeoi is used for both amp and phas
  end
  phasts = ft_specest_mtmconvol(data.trial{trialsel2},data.time{trialsel2},'freqoi',cfg.phasfreq,'timwin',cfg.phasfreqtimwin,'taper', cfg.phasfreqtaper,opt{:});
  if istrue(cfg.amptstohfats)
    ampts = hfa.trial{itrial};
  else
    ampts = ft_specest_mtmconvol(data.trial{trialsel1},data.time{trialsel1},'freqoi',cfg.ampfreq ,'timwin',cfg.ampfreqtimwin ,'taper', cfg.ampfreqtaper, opt{:});
    ampts = abs(ampts);
  end
  
  % failsafe
  if strcmp(cfg.handlespecedges,'conservative')
    if any(isnan(ampts(:))) || any(isnan(phasts(:)))
      error('timeoi determination for conservative handling of time points was not successful')
    end
  end
  
  % loop over freqs and calc wPLFs
  for iampfreq = 1:nampfreq
    for iphasfreq = 1:nphasfreq
      
      % select amplitude frequency and phase frequency time series
      currampts  = ampts(:,:,iampfreq);
      currphasts = phasts(:,:,iphasfreq);
      
      % remove NaNs if liberal (there should be none otherwise)
      if strcmp(cfg.handlespecedges,'liberal')
        nonnans    = ~logical(isnan(currampts(1,:)) + isnan(currphasts(1,:)));
        currampts  = currampts(:,nonnans);
        currphasts = currphasts(:,nonnans);
      end
      
      % set n and save dof
      ntime = size(currampts,2);
      dof(itrial,iampfreq,iphasfreq) = ntime;
      
      % failsafe
      if ntime==0
        error(['no time points without NaNs for phase @' num2str(cfg.phasfreq(iphasfreq)) 'Hz and amplitude @' num2str(cfg.ampfreq(iampfreq)) 'Hz'])
      end
        
      % center ampts if requested
      if istrue(cfg.centerampfreq)
        currampts = currampts - repmat(mean(currampts,2),[1 ntime]);
      end
      
      % remove amplitudes from phasts if requested
      if istrue(cfg.phasfreqnoamp)
        currphasts = currphasts ./ abs(currphasts);
      end
      
      %%% Surrogate by circularly shifting %%%
      if dosurrbytime
        nsteps  = round(rand(nchan,2) .* (ntime-1))+1;
        % each channel needs to shifted separately
        for ichan = 1:nchan
          currampts(ichan,:)  = circshift(currampts(ichan,:), nsteps(ichan,1),2);
          currphasts(ichan,:) = circshift(currphasts(ichan,:),nsteps(ichan,2),2);
        end
      end
      %%% Surrogate by circularly shifting %%%
      
      % calc wPLFs
      currwplfs =  currampts * currphasts.';
      
      % perform normalization specific actions
      switch cfg.normalization
        case 'overtrials'
          % save sums of squares to norm wPLFs after computing inner products of all trials
          ssqampts(:,iampfreq,iphasfreq)  = ssqampts(:,iampfreq,iphasfreq)  + sum(abs(currampts).^2,2);
          ssqphasts(:,iampfreq,iphasfreq) = ssqphasts(:,iampfreq,iphasfreq) + sum(abs(currphasts).^2,2);
        case 'withintrials'
          % compute sums of squares to norm wPLFs now
          ssqampts  = sum(abs(currampts).^2,2);
          ssqphasts = sum(abs(currphasts).^2,2);
          currwplfs = currwplfs ./ (sqrt(ssqampts) * sqrt(ssqphasts).');
          % if keeptrials = no, divide wplfs by ntrial such that the below sum over trials works for both types of normalization
          if ~istrue(cfg.keeptrials)
            currwplfs = currwplfs ./ ntrial;
          end
      end
      
      % insert into wplf
      if ~istrue(cfg.keeptrials)
        wplf(:,:,iampfreq,iphasfreq) = wplf(:,:,iampfreq,iphasfreq) + currwplfs;
      else
        wplf(itrial,:,:,iampfreq,iphasfreq) = currtrialwplfs; 
      end
    end
  end
end


% normalize wPLFs if normalization is overtrials
switch cfg.normalization
  case 'overtrials'
    if ~istrue(cfg.keeptrials)
      for iampfreq = 1:nampfreq
        for iphasfreq = 1:nphasfreq
          % construct current normalization factor
          normfactor = sqrt(ssqampts(:,iampfreq,iphasfreq)) * sqrt(ssqphasts(:,iampfreq,iphasfreq)).';
          % norm
          wplf(:,:,iampfreq,iphasfreq) = wplf(:,:,iampfreq,iphasfreq) ./ normfactor;
        end
      end
    else
      for itrial = 1:ntrial
        % extract current trial
        currwplf = squeeze(wplf(itrial,:,:,iampfreq,iphasfreq));
        for iampfreq = 1:nampfreq
          for iphasfreq = 1:nphasfreq
            % construct current normalization factor
            normfactor = sqrt(ssqampts(:,iampfreq,iphasfreq)) * sqrt(ssqphasts(:,iampfreq,iphasfreq)).';
            % norm
            currwplf(:,:,iampfreq,iphasfreq) = currwplf(:,:,iampfreq,iphasfreq) ./ normfactor;
          end
        end
        % insert normalized wplfs
        wplf(itrial,:,:,iampfreq,iphasfreq) = currwplf;
      end
    end
  case 'withintrials'
    % do nothing
end


































