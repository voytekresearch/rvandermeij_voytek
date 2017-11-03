function fourier = rmr_SPACEinput_electrophysiology(data,freqoi,timeoi,timwin,taper,tapsmofrq)

% RMR_SPACEINPUT_ELECTROPHYSIOLOGY computes a 4-way channel-by-frequency-by-epoch-by-taper
% array of Fourier coefficients from electrophysiological potential time-series to be used
% for SPACE decompositions. 
%
%  Use as:
%  fourier = rmr_SPACEinput_electrophysiology(data,freqoi,timeoi,timwin,taper,tapsmofrq)
%
% This function convolves electrophysiologicalpotential time-series with wavelets to create 
% Fourier coefficients. These Fourier coefficients can be used in SPACE decompositions to find 
% networks/components described by between-electrode phase-coupling over multiple frequencies.
% After obtaining channel-by-taper (time) matrices of Fourier coefficients, CSDs are calculated, 
% and per frequency and epoch an Eigenvalue decomposition is done. Eigenvalues less than 
% tolerance are removed, and the remaining Eigenvectors are each multiplied with the square root
% of their Eigenvalue. The resulting matrices have the same cross-product as the original 
% channel-by-taper matrices of Fourier coefficients, but is typically much smaller. 
%
% This array can be used as memory-efficient input for SPACE decompositions, 
% by taking advantage of sufficiency of cross-products, which SPACE inherited from 
% PARAFAC2 (see Kiers et al 1999, JChemoMetr.).
%
% Epoch in the above can refer to any kind of temporal segmentation (such as trials), and should be 
% a temporal division for which it makes sense to obtain separate cross spectra, and, as such, a seperate
% weight in the final epoch profile of the SPACE decomposition. Importantly, the epoch-specific sections of 
% the Fourier array are normalized by epoch length in seconds. 
% 
% Tapers in the above, stands for post-convolution time-windows centered at each time point of the input data. 
% These time-points are called tapers here, because a more complicated scheme than cross-products over time can be 
% used to compute the Fourier array as input for SPACE decompositions, such as Slepian/DPSS tapering. 
% 
% Because the relevant computions in SPACE-time/FSP are always performed on double 
% precision numbers, the output of this function is SINGLE precision. 
%
% Complex wavelet convolution is performed by multiplication in the frequency domain. Each wavelet is a 
% complex exponential (complex sinusoid) of a certain length tapered by a window function. Length is advised to 
% be dependent on cycle length of each frequency of interest.
%
%
%       Input:
%            data = FieldTrip raw data structure as obtained from the FT_PREPROCESSING function
%          freqoi = 1xNfreq vector, center frequencies of interest to use for wavelet convolution
%          timwin = 1xNfreq vector, wavelet length in seconds
%          timeoi = 1xNtime vector, time-points of interest, i.e. samples of the convolved time-series
%                                   to use for computing cross spectra.
%                                   Note: these time-points should be non-NaN for every frequency
%           taper = string, 'hanning' or 'dpss', tapering function to use in creating wavelet
%       tapsmofrq = optional, 1xNfreq vector, frequency-specific half-width smoothing box when using DPSS tapering
%
%      Output:
%         fourier = 4-way chan_freq_epoch_taper array of fourier coefficients (SINGLE precision)
%
%
%
% Copyright (C) 2016-present, Roemer van der Meij


% set n's
ntrial    = numel(data.trial);
nchan     = numel(data.label);
nfreq     = numel(freqoi);

% do some sanity checks
if nfreq~=numel(timwin)
  error('improper size of freqoi or timwin')
end
if strcmp(taper,'dpss') && (~exist('tapsmofrq','var') || isempty(tapsmofrq) || numel(tapsmofrq)~=nfreq)
  error('need to define proper tapsmofrq when using dpss tapering')
end

% get ntaper estimate
switch taper
  case 'hanning'
    ntaper = min([nchan numel(timeoi)]);
  case 'dpss'
    ntaper = min([nchan numel(timeoi)*max((2*timwin.*tapsmofrq)-1)]);
  otherwise
    error('taper not supported')
end

% pre-allocate
fourier = complex(NaN(nchan,nfreq,ntrial,ntaper,'single'),NaN(nchan,nfreq,ntrial,ntaper,'single'));
padding = ceil(max(cellfun(@numel,data.time)) ./ data.fsample);
if mod(padding,2)~=0
  padding = padding + 1; % to ensure a half hz resolution
end

% get fourier coefficients using specest
for itrial = 1:ntrial
  
  % call spectest, loop over freq
  disp(['getting fourier coefficients of trial #' num2str(itrial)])
  currtrialdat  = data.trial{itrial};
  currtrialtime = data.time{itrial};
  for ifreq = 1:nfreq
    % switch over tapers
    switch taper
      
      case 'hanning'
        tapopt = {'taper','hanning'};
        
      case 'dpss'
        % use slepian tapering only when corrected shannon number is >1 (the first 1 is 1 sec taper)
        if ((2*timwin(ifreq)*tapsmofrq(ifreq))-1) > 1
          tapopt = {'taper','dpss','tapsmofrq',tapsmofrq(ifreq)};
        else
          warning('number of Slepian tapers is 1, using Hanning taper instead')
          tapopt = {'taper','hanning'};
        end
        
      otherwise
        error('taper not supported')
    end
    
    % get spec-est
    [spectrum,ntaperout,freqdum,timeoiout] = ft_specest_mtmconvol(currtrialdat, currtrialtime,'freqoi',freqoi(ifreq),'timeoi',timeoi,'timwin',timwin(ifreq),'pad',padding,tapopt{:},'verbose',0,'polyorder',1);
    
    % reorganize to into a chan-by-tap*time matrix
    spectrum = permute(spectrum,[2 1 3 4]);
    spectrum = spectrum(:,:); % this unfolds the dimensions other than chan

    % remove NaNs (always the same over channels)
    spectrum(:,isnan(spectrum(1,:))) = [];
    
    % undo specest_mtmconvol taper-length normalization
    spectrum = spectrum ./ sqrt(2 ./ round(timwin(ifreq).*data.fsample));
    
    % correct for effective number of tapers and time-points
    spectrum = spectrum ./ sqrt(ntaperout.*numel(timeoiout)); % number of tapers is the same for each time point
    
    % compute the csd over tapers and over windows
    csd = spectrum*spectrum';
    
    % CSD failsafes
    if any(isnan(csd(:)))
      error(['NaNs remaining at ' num2str(freqoi(ifreq)) 'Hz with wavelet length of ' num2str(timwin(ifreq)) 's for trial ' num2str(itrial)])
    end
    if isempty(csd) || all(csd(:)==0)
      error(['spectral analysis at ' num2str(freqoi(ifreq)) 'Hz with wavelet length of ' num2str(timwin(ifreq)) 's resulted in zero non-NaN time points for trial ' num2str(itrial)])
    end
    
    % Reduce SPACE memory load and computation time by replacing each chan_taper matrix by the
    % Eigenvectors of its chan_chan cross-products weighted by sqrt(Eigenvalues).
    % This possible because (1) SPACE only uses the cross-products of the chan_taper matrices
    % (i.e. the frequency- and trial-specific CSD) and (2) the Eigendecomposition of a symmetric 
    % matrix A is A = VLV'. 
    % As such, VL^.5 has the same cross-products as the original chan_taper matrix.
    [V L] = eig(csd);
    L     = diag(L);
    tol   = max(size(csd))*eps(max(L)); % compute tol using matlabs default
    zeroL = L<tol;
    eigweigth = V(:,~zeroL)*diag(sqrt(L(~zeroL)));
    
    % positive semi-definite failsafe
    if any(L<-tol)
      error(['some Eigenvalues were negative: min eig: ' num2str(L(find(L==min(L),1))) ', max eig: ' num2str(L(find(L==max(L),1))) ', tol: ' num2str(tol) ' at ' num2str(freqoi(ifreq)) 'Hz and trial ' num2str(itrial)])
    end
    if  any(~isreal(L))
      error(['some Eigenvalues were complex at ' num2str(freqoi(ifreq)) 'Hz and trial ' num2str(itrial)])
    end 
      
    % save in fourier
    currm = size(eigweigth,2);
    fourier(:,ifreq,itrial,1:currm) = eigweigth;
  end
  clear spectrum
end

% safety check, just in case
if ~all([nchan,nfreq,ntrial,ntaper] == size(fourier))
  error('something went wrong, fourier has grown')
end

% trim fourier (actual ntaper can be lower than ntaper, which depends on the data, hence the need for trimming)
notnan = logical(squeeze(sum(sum(~isnan(squeeze(fourier(1,:,:,:))),2),1)));
fourier = fourier(:,:,:,notnan);














