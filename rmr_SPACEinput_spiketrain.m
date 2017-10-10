function fourier = rmr_SPACEinput_spiketrain(dat,fsample,freqoi,timwin,method)

% RMR_SPACEINPUT_SPIKETRAIN computes a 4-way channel-by-frequency-by-epoch-by-taper
% array of Fourier coefficients from neural spiking time series to be used
% for SPACE decompositions. 
%
%  Use as:
%  fourier = rmr_SPACEinput_spiketrain(dat,freqoi,timwin,taper,tapsmofrq)
%
% This function convolves binary spike trains with complex exponentials (non-tapered wavelets)
% to create 'fake' Fourier coefficients. These Fourier coefficients can be used in a SPACE-time
% decomposition to find networks of neurons defined by consistent spike timing between them.
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
% weight in the final epoch profile of the SPACE decomposition.  Importantly, the epoch-specific sections of 
% the Fourier array are normalized by epoch length in seconds. 
% 
% Tapers in the above, stands for post-convolution time-windows centered at each time point of the input data. 
% These time-points are called tapers here, because a more complicated scheme than cross-products over time can be 
% used to compute the Fourier array as input for SPACE decompositions, such as Slepian/DPSS tapering. This is not
% relevant for a Fourier array computed from spike trains however.
% 
% Because the relevant computions in SPACE-time/FSP are always performed on double 
% precision numbers, the output of this function is SINGLE precision. 
%
% The complex exponentials (complex sinusoids) used for convolution are not tapered by a window function. This is
% necessary to optimize sensitivity of the cross spectra. Their length should reflect the (broad) range at which 
% spike timing consistency is expected to be present. 
% 
%
%       Input:
%             dat = 1xNepoch cell-arrays of binary Nneuron-by-Nsample spike trains (can be sparse)
%         fsample = scalar, sampling rate in Hz
%          freqoi = 1xNfreq vector, frequencies at which to perform complex exponential convolution
%          timwin = scalar, length in seconds of the complex exponentials
%          method = 'sparseconv', 'segmentsparseconv' or 'segmentfft', method to use for convolving spike trains 
%                       - sparseconv: use convolution of sparse matrices in their original size
%                       - segmentsparseconv: same as sparseconv, but less memory intensive 
%                                             (cuts epochs into segments, taking care of overlap)
%                       - segmentfft: convolution via multiplication in the frequency domain, more memory
%                                             efficient when using a long time window (e.g. >100ms)
%
%       Output:
%         fourier = 4-way channel-by-frequency-by-epoch-by-taper array of Fourier coefficients (SINGLE precision)
%
%
%
% Copyright (C) 2016-present, Roemer van der Meij
%

% check input
if nargin == 4
  method = 'segmentsparseconv';
end
if any(strcmp(method,{'sparseconv','segmentsparseconv'}))
  if ~exist('sconv2','file')
    error('sconv2.m required from MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/41310')
  end
end
if ~iscell(dat)
  error('dat needs to a cell-array of binary spiking time series')
end

% set n's
nepoch    = numel(dat);
nchan     = size(dat{1},1);
nfreq     = numel(freqoi);

% do some sanity checks
if numel(timwin)~=1
  error('improper size of timwin')
end
if any(max(cellfun(@max,cellfun(@size,dat,'uniformoutput',0))) < round(timwin*fsample))
  error('length of complex exponential cannot exceed minimum epoch length')
end

% pre-allocate
fourier = complex(NaN(nchan,nfreq,nepoch,nchan,'single'),NaN(nchan,nfreq,nepoch,nchan,'single'));

% get fourier coefficients using specest
for iepoch = 1:nepoch
  
  % convolve with wavelets per freq
  currepochdat  = double(sparse(dat{iepoch})); % ensure sparse and double
  nsmpcurrepoch = size(currepochdat,2);
  for ifreq = 1:nfreq
    disp(['getting fourier coefficients of epoch #' num2str(iepoch) ' @' num2str(freqoi(ifreq)) 'Hz'])
    
    % construct wavelet and convolve using sparse convolution (from the file exchange)
    nsmptimwin = round(timwin .* fsample);
    anglein    = (-(nsmptimwin-1)/2 : (nsmptimwin-1)/2)   .*  ((2.*pi./fsample) .* freqoi(ifreq));
    wlt        = complex(cos(anglein), sin(anglein));
   
    %     %% debug plotting
    %         if iepoch == 1
    %           figure('name',['taper #' num2str(ifreq) ' @ ' num2str(freqoi(ifreq)) 'Hz' ],'NumberTitle','off');
    %           plot(real(wlt));
    %           hold on;
    %           plot(imag(wlt),'color','r');
    %           plot(abs(wlt),'color','gr');
    %           %plot(angle(wlt) ./ ((2*pi) ./ max(abs(wlt))),'color',[.5 .5 .5]);
    %           legend('real','imag','abs');
    %         end
    %     %% debug plotting
    
    % switch methods
    switch method
      
      
      case {'segmentfft','segmentsparseconv'}
        
        % Convolve data with wavelet and compute CSD
        % Though the spike train is sparse, the resulting spectrum is not. As such, it becomes very large very quickly
        % Therefore, calculate the csd as a running sum, running over time-windows
        nsmpwin = 5e6; % should require a couple of GB with a couple of hundred units
        nsmpwin = max([nsmpwin nsmptimwin]); % ensure window is at least as long as wavelet
        % create indices for windows, overlapping at least half of the wavelet length (if spikes occur at end of
        allind = 1:nsmpwin:nsmpcurrepoch;
        if allind(end) ~= nsmpcurrepoch
          if numel(allind)==1 % epoch is shorter than nsmpwin
            allind(end+1) = nsmpcurrepoch+1; % the +1 precorrects for the subtraction below
          elseif (nsmpcurrepoch-allind(end)) > ceil(nsmptimwin*1.5)
            % difference to end of last window is of sufficient size
            allind(end+1) = nsmpcurrepoch+1; % the +1 precorrects for the subtraction below
          else
            % difference to end of last window is of insufficient size, change to endsample
            allind(end) = nsmpcurrepoch+1;  % the +1 precorrects for the subtraction below
          end
        end
        winind = [];
        winind(:,1) = allind(1:end-1);
        winind(:,2) = allind(2:end)-1;
        nsmpoverl = ceil(nsmptimwin/2);
        nwin = size(winind,1);
        
        % loop over windows and calculate csd as running sum
        csd = zeros(nchan,nchan);
        for iwin = 1:nwin
          begind = winind(iwin,1);
          endind = winind(iwin,2);
          
          % add overlap
          if nwin == 1
            % do nothing
          elseif iwin == 1
            endind = endind + nsmpoverl;
          elseif iwin == nwin
            begind = begind - nsmpoverl;
          else
            begind = begind - nsmpoverl;
            endind = endind + nsmpoverl;
          end
          
          % calc spectrum
          switch method
            case 'segmentfft'
              
              % linear convolution using fft
              nsmpcurrwin   = numel(begind:endind);
              spctrmcurrwin = fft(cat(2,full(currepochdat(:,begind:endind)),zeros(nchan,nsmptimwin)), [], 2);
              spctrmwlt     = fft([wlt zeros(1,nsmpcurrwin)]);
              spectrum      = ifft(bsxfun(@times,spctrmcurrwin,spctrmwlt), [], 2);
              clear spctrmwlt spctrmcurrwin
              spectrum      = spectrum(:,floor(nsmptimwin/2)+1:end-ceil(nsmptimwin/2)); % remove padding
              
            case 'segmentsparseconv'
              
              % get spectrum via sparse convolution
              spectrum = sconv2(currepochdat(:,begind:endind),wlt,'same');
              
          end
          
          % remove overlap
          if nwin == 1
            ind = 1 : numel(begind:endind);
          elseif iwin == 1
            ind = 1 : (numel(begind:endind) - nsmpoverl);
          elseif iwin == nwin
            ind = (nsmpoverl + 1) : numel(begind:endind);
          else
            ind = (nsmpoverl + 1) : (numel(begind:endind) - nsmpoverl);
          end
          spectrum = spectrum(:,ind);
          
          % add a tiny bit of noise to the csd, at channel level, in case one of the channels has 0 spikes (better to add uniformily, i.e. over all channels)
          spectrum(:,1) = spectrum(:,1) +  (wlt(1).*exp(1i*2*pi*rand(nchan,1)).*eps);
          
          % calc csd as running sum
          csd = csd + spectrum*spectrum';
          clear spectrum
        end
        
        
      case 'sparseconv'
        
        % get spectrum via sparse convolution
        spectrum = sconv2(currepochdat,wlt,'same');
        
        % add a tiny bit of noise to the csd, at channel level, in case one of the channels has 0 spikes (better to add uniformily, i.e. over all channels)
        spectrum(:,1) = spectrum(:,1) +  (wlt(1).*exp(1i*2*pi*rand(nchan,1)));
        
        % calc csd
        csd = full(spectrum*spectrum'); % ensure full
        
        
      otherwise
        error('method not supported')
        
    end
    
    
    % correct csd for number of time-points (for variable length epochs)
    csd = csd ./ (size(currepochdat,2) ./ fsample); % convert to coupling/second
    
    % NaN failsafe
    if any(isnan(csd(:)))
      error(['NaNs remaining at ' num2str(freqoi(ifreq)) 'Hz and trial ' num2str(itrial)])
    end
    if isempty(csd) || all(csd(:)==0)
      error(['complex exponential convolution at ' num2str(freqoi(ifreq)) 'Hz with window length of ' num2str(timwin) 's resulted in zero non-NaN time points for trial ' num2str(itrial)])
    end
    
    
    % Reduce SPACE memory load and computation time by replacing each chan_taper matrix by the
    % Eigenvectors of its chan_chan cross-products weighted by sqrt(Eigenvalues).
    % This possible because (1) SPACE only uses the cross-products of the chan_taper matrices
    % (i.e. the frequency- and epoch-specific CSD) and (2) the Eigendecomposition of a symmetric
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
    fourier(:,ifreq,iepoch,1:currm) = eigweigth;
  end
  clear spectrum
end

% safety check, just in case
if ~all([nchan,nfreq,nepoch,nchan] == size(fourier))
  error('something went wrong, fourier has grown')
end

% trim fourier (actual ntaper can be lower than ntaper, which depends on the data, hence the need for trimming)
notnan = logical(squeeze(sum(sum(~isnan(squeeze(fourier(1,:,:,:))),2),1)));
fourier = fourier(:,:,:,notnan);









