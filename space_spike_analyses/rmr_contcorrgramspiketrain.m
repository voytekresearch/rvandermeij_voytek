function [contcorrgram,timeout] = rmr_contcorrgramspiketrain(cfg,data)

% RMR_CONTCORRGRAMSPIKETRAIN
%
%  Use as:
%  data = rmr_corrgramspiketrain(cfg,data)
%
% This function computes between-unit "correlograms" using spike trains as input.
% This is not a true histogram, but a 'continuous' version, where each count is smoothed
% with a Gaussian of cfg.gaussfwhm FWHM (...spike density)
%
% Input data should be a FieldTrip-style data structure modified by Roemer to hold
% spike-train data.
%
%
%      Cfg-options:
%          cfg.range         = 2x1 vector, time range in seconds for constructing correlograms
%          cfg.keeptrials    = 'yes' or 'no'
%          cfg.gaussfwhm     = FWHM in seconds
%          cfg.fsample       = sampling rate to save corrgram at (must be integer of sampling rate)
%

% input checks
ft_checkconfig(cfg,'required',{'range'});
cfg.keeptrials  = ft_getopt(cfg, 'keeptrials',   'no');
cfg.gaussfwhm   = ft_getopt(cfg, 'gaussfwhm',     .5 / 1000);
cfg.fsample     = ft_getopt(cfg, 'fsample',       4000);
cfg.keeptrials  = istrue(cfg.keeptrials);
if rem(data.fsample,cfg.fsample)~=0
  error('cfg.fsample needs to be integer multiple of data.fsample');
end

% set
nunit    = numel(data.label);
ntrial   = numel(data.trial);
rangesmp = round(cfg.range .* data.fsample);
time     = (rangesmp(1):rangesmp(2)) ./data.fsample;
downfac  = data.fsample/cfg.fsample;
timeout  = time(1:downfac:end);
ntime    = numel(timeout);


% make gaussian used for smoothing counts
gausssig  = cfg.gaussfwhm ./ (2.*sqrt(-2.*log(.5)));
gausscoef = normpdf(time,0,gausssig);
gausscoef = gausscoef ./ max(gausscoef);
% for speedup
%gausscoef(gausscoef<=1e-6) = 0;
gausscoef(time<-(cfg.gaussfwhm*2) | time>(cfg.gaussfwhm*2)) = 0;


%%%%% First, compute spike correlograms (looping over spikes)
if cfg.keeptrials
  contcorrgram = zeros(ntrial,nunit,nunit,ntime,'single');
else
  contcorrgram = zeros(nunit,nunit,ntime,'single');
end

% loop over trials, and count
for itrial = 1:ntrial
  disp(['getting corrgram for ' num2str(itrial) '/' num2str(ntrial)])
  
  % obtain unit id's and sample id's of spikes
  currtrial = data.trial{itrial};
  nsample   = size(currtrial,2);
  [unitid, smpid] = ind2sub(size(currtrial),find(currtrial));
  
  % loop over spikes, and count co-occurances
  for ispike = 1:numel(smpid)
    
    % subselect spikes to count
    seedsmpid  = smpid(ispike);
    spikesel   = (smpid >= (seedsmpid+rangesmp(1))) & (smpid <= (seedsmpid+rangesmp(2)));
    seedunitid = unitid(ispike);
    targunitid = unitid(spikesel);
    unitarg    = unique(targunitid);
    ntargunits = numel(unitarg);
    
    % select spike train
    selind = [seedsmpid+rangesmp(1) seedsmpid+rangesmp(2)];
    if selind(1)<1
      prepad    = sparse([],[],[],ntargunits,abs(selind(1))+1);
      selind(1) = 1;
    else 
      prepad = [];
    end
    if selind(2)>nsample
      postpad   = sparse([],[],[],ntargunits,selind(2)-nsample);
      selind(2) = nsample;
    else
      postpad = [];
    end
    % pad to full length
    currtrain = currtrial(unitarg,selind(1):selind(2));
    currtrain = [prepad currtrain postpad];
    counts    = sconv2(currtrain,gausscoef,'same');
    % 'downsample'
    counts = counts(:,1:downfac:end);
    
    % count
    if cfg.keeptrials
      contcorrgram(itrial,seedunitid,unitarg,:) = permute(contcorrgram(itrial,seedunitid,unitarg,:),[3 4 1 2]) + full(counts);
    else
      contcorrgram(seedunitid,unitarg,:) = permute(contcorrgram(seedunitid,unitarg,:),[2 3 1]) + full(counts);
    end
  end
end














