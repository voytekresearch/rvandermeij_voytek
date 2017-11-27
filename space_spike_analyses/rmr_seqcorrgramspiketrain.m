function seqcorrgram = rmr_seqcorrgramspiketrain(cfg,data)

% RMR_CORRGRAMSPIKETRAIN
%
%  Use as:
%  data = rmr_corrgramspiketrain(cfg,data)
%
% This function computes between-unit correlograms using spike trains as input.
%
%
% Input data should be a FieldTrip-style data structure modified by Roemer to hold
% spike-train data.
%
%
%      Cfg-options:
%          cfg.binedges   = 1xNbins+1 vector of bin EDGES in seconds to construct correlograms
%          cfg.keeptrials = 'yes' or 'no'
%          cfg.selectunit = 1xNunit cell-array containing labels of units to search for
%
%

% input checks
ft_checkconfig(cfg,'required',{'binedges'});
cfg.keeptrials  = ft_getopt(cfg, 'keeptrials',   'no');
cfg.keeptrials  = istrue(cfg.keeptrials);

% select units
[dum, selind] = match_str(cfg.selectunit,data.label);
for itrial = 1:numel(data.trial)
  data.trial{itrial} = data.trial{itrial}(selind,:);
end
data.label = data.label(selind);

% set
nunit   = numel(cfg.selectunit);
ntrial  = numel(data.trial);
nbins   = numel(cfg.binedges)-1;

% get some vars to use later
timebinedges  = cfg.binedges;
timebincenter = cfg.binedges(2:end)-mean(diff(cfg.binedges))/2;
smpbinedges   = round(timebinedges .* data.fsample);
smpbincenter  = round(timebincenter .* data.fsample);

% prep gathering
seqcorrgram = cell(1,nunit); %1 unit spiking (empty), 2 spiking, 3 spiking, etc
%%%%% compute spike correlograms
for iunit = 2:nunit
  if cfg.keeptrials
    seqcorrgram{iunit} = zeros(ntrial,nunit,nunit,nbins);
  else
    seqcorrgram{iunit} = zeros(nunit,nunit,nbins);
  end
end

% loop over trials, and count
for itrial = 1:ntrial
  disp(['getting corrgram for ' num2str(itrial) '/' num2str(ntrial)])
  
  % obtain unit id's and sample id's of spikes
  [unitid, smpid] = ind2sub(size(data.trial{itrial}),find(data.trial{itrial}));
  
  % loop over spikes, and count co-occurances
  for ispike = 1:numel(smpid)
    
    % subselect spikes to count
    spikesel   = smpid >= (smpid(ispike)-abs(smpbinedges(1))) & smpid <= (smpid(ispike)+smpbinedges(end));
    seedsmpid  = smpid(ispike);
    seedunitid = unitid(ispike);
    targsmpid  = smpid(spikesel) - seedsmpid;
    targunitid = unitid(spikesel);
    targsmpid(targunitid==seedunitid)  = [];
    targunitid(targunitid==seedunitid) = [];
    unitarg    = unique(targunitid);
    
    % possible sequence length
    nseqact = numel(unitarg)+1;
    
    % count
    for itarg = 1:numel(unitarg)
      currtarg   = unitarg(itarg);
      currspikesel = targunitid==currtarg;
      if cfg.keeptrials
        seqcorrgram{nseqact}(itrial,seedunitid,currtarg,:) = permute(seqcorrgram{nseqact}(itrial,seedunitid,currtarg,:),[4 1 2 3])' + histcounts(targsmpid(currspikesel),smpbinedges);
      else
        seqcorrgram{nseqact}(seedunitid,currtarg,:) = permute(seqcorrgram{nseqact}(seedunitid,currtarg,:),[3 1 2])' + histcounts(targsmpid(currspikesel),smpbinedges);
      end
    end
  end
end












