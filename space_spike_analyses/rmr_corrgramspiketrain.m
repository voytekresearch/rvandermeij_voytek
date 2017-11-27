function [corrgram,timebincenter] = rmr_corrgramspiketrain(cfg,data)

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
%
%

% input checks
ft_checkconfig(cfg,'required',{'binedges'});
cfg.keeptrials  = ft_getopt(cfg, 'keeptrials',   'no');
cfg.keeptrials  = istrue(cfg.keeptrials);

% set
nunit   = numel(data.label);
ntrial  = numel(data.trial);
nbins   = numel(cfg.binedges)-1;

% get some vars to use later
timebinedges  = cfg.binedges;
timebincenter = cfg.binedges(2:end)-mean(diff(cfg.binedges))/2;
smpbinedges   = round(timebinedges .* data.fsample);
smpbincenter  = round(timebincenter .* data.fsample);

% choose between two approaches for getting the corrgram, not sure yet which is faster (overspikes is prolly faster in real situation)
method = 'overspikes'; %'overunits' 'overspikes'
switch method
  
  case 'overunits'
    %%%%% First, compute spike correlograms (looping over units)
    if cfg.keeptrials
      corrgram = zeros(ntrial,nunit,nunit,nbins);
    else
      corrgram = zeros(nunit,nunit,nbins);
    end
    
    % loop over trials, and count
    for itrial = 1:ntrial
      disp(['getting corrgram for ' num2str(itrial) '/' num2str(ntrial)])
      
      % get time stamps of all units
      spikets = cell(1,nunit);
      for iunit = 1:nunit
        spikets{iunit} = find(data.trial{itrial}(iunit,:));
      end
      
      % loop over seed units
      for iseed = 1:nunit
        currseedspk = spikets{iseed};
        
        % loop over target units
        for itarg = 1:nunit
          currtargsp = spikets{itarg};
          
          % then, from the perspective of each spike of the seed unit, count occurances of spikes of target unit in binedges
          currcount = zeros(1,nbins);
          for ispseed = 1:numel(spikets{iseed})
            spkdiffts = currtargsp - currseedspk(ispseed);
            currcount = currcount + histcounts(spkdiffts,smpbinedges);
          end
          if cfg.keeptrials
            corrgram(itrial,iseed,itarg,:) = permute(corrgram(itrial,iseed,itarg,:),[4 1 2 3]) + currcount';
          else
            corrgram(iseed,itarg,:) = permute(corrgram(iseed,itarg,:),[3 1 2]) + currcount';
          end
        end % itarg
      end % iseed
      
    end % itrial
    %%%%%
    
  case 'overspikes'
    %%%%% First, compute spike correlograms (looping over spikes)
    if cfg.keeptrials
      corrgram = zeros(ntrial,nunit,nunit,nbins);
    else
      corrgram = zeros(nunit,nunit,nbins);
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
        targsmpid  = smpid(spikesel) - seedsmpid;
        seedunitid = unitid(ispike);
        targunitid = unitid(spikesel);
        unitarg    = unique(targunitid);
        
        % count
        for itarg = 1:numel(unitarg)
          currtarg   = unitarg(itarg);
          currspikesel = targunitid==currtarg;
          if cfg.keeptrials
            corrgram(itrial,seedunitid,currtarg,:) = permute(corrgram(itrial,seedunitid,currtarg,:),[4 1 2 3])' + histcounts(targsmpid(currspikesel),smpbinedges);
          else
            corrgram(seedunitid,currtarg,:) = permute(corrgram(seedunitid,currtarg,:),[3 1 2])' + histcounts(targsmpid(currspikesel),smpbinedges);
          end
        end
      end
    end
    %%%%%
end













