function hfa = rmr_highfrequencyactivity(data,freqoi,timeoi,timwin,taper)

% RMR_HIGHFREQUENCYACTIVITY calculates High Frequency Activity time-series.
%
%  Use as:
%  hfa = rmr_highfrequencyactivity(data,freqoi,timeoi,timwin,taper)
%
% This functions calculates HFA time series by Z-scoring the amplitude of each frequency 
% in freqoi and then averaging these Z-scored amplitudes over frequencies. Spectral analysis
% is performed by convolving each trial with a complex exponential (wavelet) at each requested 
% frequency, only the center frequency of the resulting Fourier coefficients is used. The 
% complex exponential is tapered by the function given in the input. Convolution is performed by 
% multiplication in the frequency domain.
%
% To save memory, HFA time series are obtained in two passes, the first to obtain 
% the mean and std-dev of each frequency, the second to obtain the HFA time series as a running
% average. For spectral analysis, data is zero-padded out to an integer number
% of seconds.
%
%
%        Input:
%            data = FieldTrip data structure as obtained from FT_PREPROCESSING 
%          freqoi = 1xNfreq vector of frequencies to use for calculating HFA
%          timeoi = 'all', or 1xNtime vector of time-points of interest 
%          timwin = scalar, sliding time-window length in seconds to use for spectral analysis
%           taper = string, 'hanning' or other window
%
%        Output:
%             hfa = a FieldTrip data structure as obtained from FT_PREPROCESSING, but with HFA time series instead of voltages
%
%
% This function depends on the FieldTrip toolbox, which can be found at www.fieldtriptoolbox.org.
%
% Copyright (C) 2016-present, Roemer van der Meij
%

% check inputs
data = ft_checkdata(data, 'datatype', {'raw', 'raw+comp', 'mvar'}, 'hassampleinfo', 'yes');
if numel(freqoi)==1 || ~isvector(freqoi)
  error('freqoi needs to be a 1xNfreq vector')
end
if (isnumeric(timeoi) && (numel(freqoi)==1 || ~isvector(timeoi) )) || (~isnumeric(timeoi) && ~strcmp(timeoi,'all'))
  error('timeoi needs to either be a 1xNtime vector or ''all''')
end  
if numel(timwin)~=1
  error('only one spectral analysis time-window size can be specified')
end
% taper is checked by specest

% set n's
ntrial    = numel(data.trial);
nchan     = numel(data.label);
nfreq     = numel(freqoi);

% expand timwin
timwin = repmat(timwin,[1 nfreq]);

% do some sanity checks
if nfreq~=numel(timwin)
  error('improper size of freqoi or timwin')
end
if strcmp(taper,'dpss') 
  error('multitapering not supported')
end

% determine padding for spectral analysis
padding = ceil(max(cellfun(@numel,data.time)) ./ data.fsample);


%%%%%%%%%%
% first, obtain mean and std-dev for each frequency
freqsum = zeros(nchan,nfreq);
freqssq = zeros(nchan,nfreq);
ntime   = 0;
for itrial = 1:ntrial
  % call spectest, loop over freq
  disp(['calculating HFA time-series: obtaining mean and std-dev for trial ' num2str(itrial) ' of ' num2str(ntrial)])
  
  % get spec-est
  spectrum = ft_specest_mtmconvol(data.trial{itrial}, data.time{itrial},'freqoi',freqoi,'timeoi',timeoi,'timwin',timwin,'taper',taper,'pad',padding,'verbose',0,'dimord','chan_time_freqtap');

  % save sums
  freqsum = freqsum + permute(nansum(abs(spectrum)   ,2),[1 3 2]);
  freqssq = freqssq + permute(nansum(abs(spectrum).^2,2),[1 3 2]);
  ntime   = ntime + sum(~isnan(permute(spectrum(1,:,1),[2 1 3])));
end
% compute mean and std-dev
freqmean = freqsum ./ ntime;
freqstd  = sqrt((freqssq - ((freqsum.^2)./ntime))./ (ntime-1));
%%%%%%%%%%


%%%%%%%%%%
% then, obtain hfa time series
hfats = cell(1,ntrial);
for itrial = 1:ntrial
  % call spectest, loop over freq
  disp(['calculating HFA time-series: trial ' num2str(itrial) ' of ' num2str(ntrial)])
  
  % get spec-est
  spectrum = ft_specest_mtmconvol(data.trial{itrial}, data.time{itrial},'freqoi',freqoi,'timeoi',timeoi,'timwin',timwin,'taper',taper,'pad',padding,'verbose',0,'dimord','chan_time_freqtap');

  % z-score
  ntime = size(spectrum,2);
  hfats{itrial} = mean((abs(spectrum) - permute(repmat(freqmean,[1 1 ntime]),[1 3 2])) ./ permute(repmat(freqstd,[1 1 ntime]),[1 3 2]),3);
end
%%%%%%%%%%

%%%%%%%%%%
% create time
time = repmat({timeoi},size(hfats));
%%%%%%%%%%

% create output
hfa              = [];
hfa.fsample      = 1./mean(diff(timeoi));
hfa.label        = data.label;
hfa.trial        = hfats;
hfa.time         = time;
hfa.trialinfo    = data.trialinfo;
hfa.cfg          = [];
hfa.cfg.freqoi   = freqoi;
hfa.cfg.timeoi   = timeoi;
hfa.cfg.timwin   = timwin;
hfa.cfg.taper    = taper;
hfa.cfg.previous = data.cfg;







