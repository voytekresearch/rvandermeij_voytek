function rmr_crcnspfc2_getfourier




% get info
info = rmr_crcnspfc2_info;
info.datapath = '/Users/roemervandermeij/Work/Data/tmpdata/CRCNS_PFC2/';

% set name suffix
%fnnamesuffix = 'four_direct_timeavg_1segsel_timwin10ms_freq50-50-1000';
fnnamesuffix = 'four_conv_timeavg_1segsel_timwin20ms_freq50-50-1000';
fnnamesuffix = 'four_conv_timeavg_1hzfiring_1segsel_timwin20ms_freq50-50-1000';

% fake fourier settings
freqoi  = 50:50:1000;
timwin  = .020;

% read data and get fourier
segsel    = 1;
% nprune    = 20;
% pruneperc = 10;
for    isess = 2%      :numel(info.session); % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  
  % read in data
  fsample = 20000;
  data = rmr_crcnspfc2_readspike([info.datapath currsess],fsample,'mazerun',0,segsel);
  % select units with firing rate above 1hz
  trllength = cat(1,data.time{:});
  trllength = trllength(:,2) - trllength(:,1);
  trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
  trlspikes = full(cat(2,trlspikes{:}));
  spikerate = mean(bsxfun(@rdivide,trlspikes,trllength.'),2);
  selind = spikerate>1;
  for itrial = 1:numel(data.trial)
    data.trial{itrial} = data.trial{itrial}(selind,:);
  end
  data.label = data.label(selind);
    
  % set fns
  fourfn = [info.savepath currsess '_' 'fourier' '_' fnnamesuffix '.mat'];
  datafn = [info.savepath currsess '_' 'dataetc' '_' fnnamesuffix '.mat'];
  if ~exist(fourfn,'file') || ~exist(datafn,'file')
    %fourier = roe_fourierspiketrain_direct(data,freqoi,timwin,true);
    fourier = rmr_SPACEinput_spiketrain(data.trial,data.fsample,freqoi,timwin,'sparseconv');
    
    % save fourier and data
    save(fourfn,'fourier','-v7.3')
    
    % save data and stuff
    data = rmfield(data,'trial');
    save(datafn,'data','freqoi','timwin')
  end
   
   
  % get pruned data sets if desired
  prunemode = {'odd','even'}; % {'odd','even'}  {'1of4','2of4','3of4','4of4'}
  for iprune = 1:numel(prunemode)
    % set fns
    currfnsuff = ['prunemode' '_' prunemode{iprune}];
    fourprfn = [info.savepath currsess '_' 'fourier' '_' fnnamesuffix '_' currfnsuff '.mat'];
    dataprfn = [info.savepath currsess '_' 'dataetc' '_' fnnamesuffix '_' currfnsuff '.mat'];
    if ~exist(fourprfn,'file') || ~exist(dataprfn,'file')
      
      % fetch data again
      data = rmr_crcnspfc2_readspike([info.datapath currsess],fsample,'mazerun',0,segsel);
      % select units with firing rate above 5hz
      trllength = cat(1,data.time{:});
      trllength = trllength(:,2) - trllength(:,1);
      trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
      trlspikes = full(cat(2,trlspikes{:}));
      spikerate = mean(bsxfun(@rdivide,trlspikes,trllength.'),2);
      selind = spikerate>1;
      for itrial = 1:numel(data.trial)
        data.trial{itrial} = data.trial{itrial}(selind,:);
      end
      data.label = data.label(selind);
      
      % prune
      datapr = rmr_prunespiketrain(data,prunemode{iprune});
      
      % get fourier
      %fourier = roe_fourierspiketrain_direct(datapr,freqoi,timwin,true);
      fourier = roe_fourierspiketrain_conv(datapr,freqoi,timwin,'sparseconv');

      
      % save fourier and data
      save(fourprfn,'fourier','-v7.3')
      
      % save data and stuff
      datapr = rmfield(datapr,'trial');
      save(dataprfn,'datapr','freqoi','timwin')
    end
  end
%   %   if nprune>0
%   %     for iprune = 1:nprune
%   %       disp(['constructing pruned dataset ' num2str(iprune) '/' num2str(nprune)])
%   %
%   %       % set fns
%   %       currfnsuff = ['pruneset' num2str(iprune) '_' num2str(pruneperc)];
%   %       fourprfn = [info.savepath currsess '_' 'fourier' '_' fnnamesuffix '_' currfnsuff '.mat'];
%   %       dataprfn = [info.savepath currsess '_' 'dataetc' '_' fnnamesuffix '_' currfnsuff '.mat'];
%   %       if ~exist(fourprfn,'file') || ~exist(dataprfn,'file')
%   %
%   %         % fetch data again
%   %         data = rmr_crcnspfc2_readspike([info.datapath currsess],fsample,'mazerun',0,segsel);
%   %
%   %         % prune
%   %         datapr = rmr_prunespiketrain(data,pruneperc);
%   %
%   %         % get fourier
%   %         %fourier = roe_fourierspiketrain_direct(datapr,freqoi,timwin,true);
%   %         fourier = roe_fourierspiketrain_conv(datapr,freqoi,timwin,'sparseconv');
%   %
%   %         % save fourier and data
%   %         save(fourprfn,'fourier','-v7.3')
%   %
%   %         % save data and stuff
%   %         datapr = rmfield(datapr,'trial');
%   %         save(dataprfn,'datapr','freqoi','timwin')
%   %       end
%   %     end
%   %   end
end


























function playground






% get info
info = rmr_crcnspfc2_info;
%info.datapath = '/Users/roemer/Work/Data/tmpdata/CRCNS_PFC2/';

% 
for    isess = 1      :numel(info.session); % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  
  % read in data
  fsample = 20000; % oversample intentionally
  data = rmr_crcnspfc2_readspike([info.datapath currsess],fsample,'mazerun',0);
  % don't select based on spikerate, this was already done (roughly 1Hz minimum)
  
  
  %%% Compute spike correlelolelolelograms
  nunit   = numel(data.label);
  ntrial  = numel(data.trial);
  
  % first set what time bins we want to use
  timebins = (-20:1:0:1:20) ./ 1000; % in s
  nbins    = numel(timebins)-1;
  
  % prepare count storage
  corrgram = zeros(nunit,nunit,nbins);
  
  % loop over trials, and count
  for itrial = 1:ntrial
    disp(['getting corrgram for ' num2str(itrial) '/' num2str(ntrial) ' of ' currsess])
    % get indices for all units, and convert to time stamps
    spikets = [];
    for iunit = 1:nunit
      spikets{iunit} = find(data.trial{itrial}(iunit,:)) ./ fsample;
    end
    % loop over seed units
    for iseed = 1:nunit
      % then, from the perspective of each spike of the seed unit, count occurances in timebins
      for ispseed = 1:numel(spikets{iseed});
        % loop over target units
        currcount = zeros(1,nbins);
        for itarg = 1:nunit
          currtargsp = spikets{itarg};
          currtargsp = currtargsp - spikets{iseed}(ispseed);
          currcount  = currcount + histcounts(currtargsp,timebins);
        end
        corrgram(iseed,itarg,:) = currcount;
      end % itarg
    end % iseed
  end % itrial

  % save
  fn = [info.savepath currsess '_' 'corrgram' '_' '-20msto20ms1msbins' '.mat'];
  save(fn,'allcorrgram','timebins')
end













