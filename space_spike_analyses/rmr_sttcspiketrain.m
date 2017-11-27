function sttc = rmr_sttcspiketrain(cfg,data,corrgram)

% RMR_STTCSPIKETRAIN
%
%  Use as:
%  data = rmr_corrgramspiketrain(cfg,data)
%
% This function computes between-unit STTC, as given by Cutts and Eglen 2014, using the corrgram as base.
%
%
% Input data should be a FieldTrip-style data structure modified by Roemer to hold
% spike-train data.
%
%
%      Cfg-options:
%          cfg.binedges   = 1xNbins+1 vector of bin EDGES in seconds to construct correlograms
%          cfg.trials     = 1xNcond cell-array containing trial indices over which to calculate the STTCs
%
%

% input checks
ft_checkconfig(cfg,'required',{'binedges'});
ft_checkconfig(cfg,'required',{'trials'});

% set
nunit   = numel(data.label);
ntrial  = numel(data.trial);
nbins   = numel(cfg.binedges)-1;
if any(cellfun(@max,cfg.trials)>ntrial)
  error('trial max woops')
end
if any(cellfun(@min,cfg.trials)<1)
  error('trial min woops')
end

% get some vars to use later
timebinedges  = cfg.binedges;
timebincenter = cfg.binedges(2:end)-mean(diff(cfg.binedges))/2;
smpbinedges   = round(timebinedges .* data.fsample);
smpbincenter  = round(timebincenter .* data.fsample);
smpperbin     = round(mean(diff(smpbinedges)));

%%%%%
sttc = zeros(numel(cfg.trials),nunit,nunit,nbins);
% per cond, for each unit pair, compute the STTC
for icond = 1:numel(cfg.trials)
  timeprop   = zeros(nunit, nbins);
  currtrials = cfg.trials{icond};
  stall      = cat(2,data.trial{currtrials});
  
  % first, compute the total spike numbers of each unit
  spikecount = full(sum(stall,2));
  
  % then, compute, for each unit and each bin, the proportion of the total recording that is encompassed by each bin (without duplication)
  nsample    = size(stall,2);
  for iunit = 1:nunit
    currst = stall(iunit,:);
    smpid  = find(currst);
    nspike = numel(smpid);
    if nspike>0
      for ibin = 1:nbins
        binsmp = cell(1,nspike);
        for ispike = 1:nspike
          ind1 = max(smpid(ispike)-round(smpperbin/2),1);
          ind2 = min(smpid(ispike)+round(smpperbin/2),nsample);
          binsmp{ispike} = ind1:ind2;
        end
        % throw out doubles
        binsmp = unique(cat(2,binsmp{:}));
        timeprop(iunit,ibin) = numel(binsmp)./nsample;
      end
    end
  end
  
  % finally, compute STTCs for each unit pair, each bin
  currcg = squeeze(sum(corrgram(currtrials,:,:,:),1));
  for iseed = 1:nunit
    for itarg = 1:nunit
      seedtp = timeprop(iseed,:).';
      targtp = timeprop(itarg,:).';
      seedcgtarg = squeeze(currcg(itarg,iseed,:)) ./ spikecount(iseed);
      targcgseed = squeeze(currcg(iseed,itarg,:)) ./ spikecount(itarg);
      currsttc = (((seedcgtarg-targtp) ./ (1 - (seedcgtarg.*targtp))) + ((targcgseed-seedtp) ./ (1 - (targcgseed.*seedtp)))) ./ 2;
      % save it
      sttc(icond,iseed,itarg,:) = currsttc;
    end
  end
end
%%%%%










