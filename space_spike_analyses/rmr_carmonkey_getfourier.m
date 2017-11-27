function rmr_carmonkey_getfourier



% get info
info = rmr_carmonkey_info;
%info.datapath = '/Users/roemervandermeij/Work/Data/tmpdata/carmena_monkey/';

% set name suffix
%fnnamesuffix = 'four_direct_timeavg_1hzfiring_timwin10ms_freq50-50-1000_gotilltar';
%fnnamesuffix = 'four_direct_timeavg_1hzfiring_timwin10ms_freq50-50-1000_centertillgo';
fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin20ms_freq50-50-1000_gotilltar';
%fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin20ms_freq50-50-1000_centertillgo';
%fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin100ms_freq10-10-200_gotilltar';


% fake fourier settings
freqoi  = 50:50:1000;
timwin  = .020;
% freqoi  = 10:10:200;
% timwin  = .100;

% extract trialdeftype
ind = strfind(fnnamesuffix,'_');
trialdeftype = fnnamesuffix(ind(end)+1:end);

% read data and get fourier
% nprune    = 20;
% pruneperc = 10;
for      isubj = 1     :numel(info.subj)
  for    isess = 1     :numel(info.session.(info.subj{isubj}))
    
    % set currs
    currsubj = info.subj{isubj};
    currsess = info.session.(currsubj){isess};
    
    % fetch data
    fsample = 40000; %
    data = rmr_carmonkey_readspike([info.datapath currsess],fsample,trialdeftype);
    
    % select units with firing rate above 1 hz
    trllength = cat(1,data.time{:});
    trllength = trllength(:,2) - trllength(:,1);
    trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
    trlspikes = full(cat(2,trlspikes{:}));
    spikerate = sum(trlspikes,2) ./ sum(trllength);
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
      fourier = roe_fourierspiketrain_conv(data,freqoi,timwin,'segmentfft'); % segmentsparseconv segmentfft sparseconv
      
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
        fsample = 40000; %
        data = rmr_carmonkey_readspike([info.datapath currsess],fsample,trialdeftype);
        
        % select units with firing rate above 1 hz
        trllength = cat(1,data.time{:});
        trllength = trllength(:,2) - trllength(:,1);
        trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
        trlspikes = full(cat(2,trlspikes{:}));
        spikerate = sum(trlspikes,2) ./ sum(trllength);
        selind = spikerate>1;
        for itrial = 1:numel(data.trial)
          data.trial{itrial} = data.trial{itrial}(selind,:);
        end
        data.label = data.label(selind);
        
        % prune
        datapr = rmr_prunespiketrain(data,prunemode{iprune});
        
        % get fourier
        %fourier = roe_fourierspiketrain_direct(datapr,freqoi,timwin,true);
        fourier = roe_fourierspiketrain_conv(datapr,freqoi,timwin,'sparseconv'); % segmentsparseconv segmentfft sparseconv
        
        
        % save fourier and data
        save(fourprfn,'fourier','-v7.3')
        
        % save data and stuff
        datapr = rmfield(datapr,'trial');
        save(dataprfn,'datapr','freqoi','timwin')
        clear fourier data datapr
      end
    end
    %   if nprune>0
    %     for iprune = 1:nprune
    %       disp(['constructing pruned dataset ' num2str(iprune) '/' num2str(nprune)])
    %
    %       % set fns
    %       currfnsuff = ['pruneset' num2str(iprune) '_' num2str(pruneperc)];
    %       fourprfn = [info.savepath currsess '_' 'fourier' '_' fnnamesuffix '_' currfnsuff '.mat'];
    %       dataprfn = [info.savepath currsess '_' 'dataetc' '_' fnnamesuffix '_' currfnsuff '.mat'];
    %       if ~exist(fourprfn,'file') || ~exist(dataprfn,'file')
    %
    %         % fetch data again
    %         data = rmr_crcnspfc2_readspike([info.datapath currsess],fsample,'mazerun',0,segsel);
    %
    %         % prune
    %         datapr = rmr_prunespiketrain(data,pruneperc);
    %
    %         % get fourier
    %         %fourier = roe_fourierspiketrain_direct(datapr,freqoi,timwin,true);
    %         fourier = roe_fourierspiketrain_conv(datapr,freqoi,timwin,'sparseconv'); % segmentsparseconv segmentfft sparseconv
    %
    %         % save fourier and data
    %         save(fourprfn,'fourier','-v7.3')
    %
    %         % save data and stuff
    %         datapr = rmfield(datapr,'trial');
    %         save(dataprfn,'datapr','freqoi','timwin')
    %       end
    %     end
    %   end
  end
end



