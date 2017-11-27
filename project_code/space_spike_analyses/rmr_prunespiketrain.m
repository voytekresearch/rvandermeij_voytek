function pruneddata = rmr_prunespiketrain(data,prunemode)

% RMR_PRUNESPIKETRAIN
%
%  Use as:
%  data = rmr_prunespiketrain(data, prunemode)
%
%  Where prunemode is a either percentage of spikes to remove from each unit at random,
%  or 'odd' or 'even'
%
% Input should be a FieldTrip style data structure (Note, not a spike data structure),
% where thee spike times are represented as a boolean spike train. Each unit is
% a separate channel.
%
% The data is segmented into trials.
% Two exceptions w.r.t. regular fieldtrip data structure.
% 1) data.trial is sparse
% 2) data.time only contains [begtime endtime]
%

% check input
if isnumeric(prunemode) && prunemode>100
  error('prunemode should be a % when specified numerically')
end
if ~isnumeric(prunemode) && ~any(strcmp(prunemode,{'odd','even','1of4','2of4','3of4','4of4'}))
    error('incorrect prunemode')
end
  


% set random seed using clock
randseed = sum(clock.*1e6);
rng(randseed);

% set
nunit   = numel(data.label);
ntrial  = numel(data.trial);
nsmptrl = cellfun(@size,data.trial,repmat({2},[1 ntrial]));

% First, concatenate spike trains and convert to indices per unit (!)
spiketrains = cat(2,data.trial{:});
smpid  = cell(1,nunit);
for iunit = 1:nunit
  smpid{iunit} = find(spiketrains(iunit,:));
end

% then, get the indices of spikes to prune
prsmpid = cell(1,nunit);
for iunit = 1:nunit
  currspikes = smpid{iunit};
  nspikes    = numel(currspikes);
  if isnumeric(prunemode)
    remind     = randperm(nspikes, round(nspikes*(prunemode./100)));
  else
    switch prunemode
      
      case 'odd'
        remind = 2:2:nspikes;
      case 'even'
        remind = 1:2:nspikes;
        
      case '1of4'
        remind = 1:4:nspikes;
      case '2of4'
        remind = 2:4:nspikes;
      case '3of4'
        remind = 3:4:nspikes;
      case '4of4'
        remind = 4:4:nspikes;
        
      otherwise
        error('unsupported prunemode')
    end
  end
  prsmpid{iunit} = sort(currspikes(remind));
end

% now, remove spikes from spike trains
prtrial = data.trial;
trlbegind = [0 cumsum(nsmptrl(1:end-1))] + 1;
trlendind = cumsum(nsmptrl);
for iunit = 1:nunit
  for itrial = 1:ntrial
    currprsmpid = prsmpid{iunit}(prsmpid{iunit}>=trlbegind(itrial) & prsmpid{iunit}<=trlendind(itrial));
    currprsmpid = currprsmpid-trlbegind(itrial) + 1;
    prtrial{itrial}(iunit,currprsmpid) = 0;
  end 
end


% create output data
pruneddata = data;
pruneddata.trial         = prtrial;
pruneddata.cfg           = [];
pruneddata.cfg.previous  = data.cfg;
pruneddata.cfg.pruneperc = prunemode;
pruneddata.cfg.randseed  = randseed;





















