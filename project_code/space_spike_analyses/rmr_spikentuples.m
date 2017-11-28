function spikentuples = rmr_spikentuples(cfg,data,nwaycomp)

% RMR_SPIKENTUPLES
%
%
% This function searches for spike sequences/ntuples, within the spike train,
% as specifief by each spiking network in nwaycomp, using its time-delay map to
% set the sequence.
%
%
%
%      Cfg-options:
%          cfg.nunit     = scalar, number of (strength-sorted) units of network to search for
%          cfg.jitter    = scalar, half-bin in seconds to search within (rounded up in samples)
%
%
%      Output: 
%           spikentuples = ncompX1 cell-array
%            {each cell} = ntupleXnmemberXnseedXntrial IJKL array of counts
%                            *ntupleXnmember matrices*  = mat(i,j==k)  counts of i-tuple at j==k seed in l trial
%                                                       = mat(i,j~=k)  counts of unit j~=k present in i-tuple at k seed in l trial
%
%


% input checks
ft_checkconfig(cfg,'required',{'nunit','jitter'});


% set
nunit     = numel(data.label);
nmember   = cfg.nunit;
ntrial    = numel(data.trial);
jittersmp = round(cfg.jitter .* data.fsample);
ncomp     = numel(nwaycomp.comp);


% parse nwaycomp
A = zeros(nunit,ncomp);
S = zeros(nunit,ncomp);
for icomp = 1:ncomp
  A(:,icomp) = nwaycomp.comp{icomp}{1};
  S(:,icomp) = nwaycomp.comp{icomp}{4};
end


% alloc
spikentuples = cell(ncomp,1);

% loop over comp and find
for icomp = 1:ncomp
  
  % get comp specifics
  [dum sortind] = sort(A(:,icomp),'descend');
  unitind    = sortind(1:nmember);
  timeseq    = S(unitind,icomp);
  seqdiffmat = [];
  for imember = 1:nmember
    seqdiffmat(imember,:) = timeseq - timeseq(imember);
  end
  seqdiffmat = round(seqdiffmat * data.fsample);
  
  % loop over trials
  trialntuples = zeros(nmember,nmember,nmember,ntrial);
  for itrial = 1:ntrial
    disp(['finding ntuples for comp ' num2str(icomp) '/' num2str(ncomp) ' in trial ' num2str(itrial) '/' num2str(ntrial)])
    
    % obtain unit id's and sample id's of spikes
    currtrial = data.trial{itrial}(unitind,:); % unitid now reflects order of strongest in network
    nsample   = size(currtrial,2);
    smpid     = cell(nmember,1);
    for imember = 1:nmember
      smpid{imember} = find(currtrial(imember,:));
    end
    
    % loop over seedunits
    unitntuples = zeros(nmember,nmember,nmember);
    for iseedu = 1:nmember
      targu = setdiff(1:nmember,iseedu);
      seedntuples = zeros(nmember);
      % loop seed spikes
      for iseeds = 1:numel(smpid{iseedu})
        
        %
        seedsmpid = smpid{iseedu}(iseeds);
        currseq   = zeros(nmember,1);
        currseq(iseedu) = 1;
        
        % search members
        for itargu = targu
          range = zeros(2,1);
          range(1) = max([1 seedsmpid+seqdiffmat(iseedu,itargu)-jittersmp]);
          range(2) = min([nsample seedsmpid+seqdiffmat(iseedu,itargu)+jittersmp]);
          if any(  (smpid{itargu} >= range(1)) & (smpid{itargu} <= range(2))  )
            currseq(itargu) = 1;
          end
        end
        
        % save
        seedntuples(sum(currseq),:) = seedntuples(sum(currseq),:) + currseq';
      end % iseeds
      unitntuples(:,:,iseedu) = seedntuples;
    end % iseedu
    trialntuples(:,:,:,itrial) = unitntuples;
  end % itrial
  spikentuples{icomp} = trialntuples;
end % icomp













