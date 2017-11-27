function rmr_crcnsssc3_getfourier


    
    
% get info
info = rmr_crcnsssc3_info;
%info.datapath = '/Users/roemervandermeij/Work/Data/tmpdata/CRCNS_SSC3/';

% set name suffix
fnnamesuffix = 'four_conv_timeavg_timwin20ms_freq50-50-1000';

% fake fourier settings
freqoi  = 50:50:1000;
timwin  = .020;

% read data and get fourier
triallength = 50;
% nprune    = 20;
% pruneperc = 10;
for    isess = 11%      :numel(info.session); % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  
  % read in data
  data = rmr_crcnsssc3_readspike([info.datapath currsess],triallength);
  
  % select units with firing rate at least 1Hz in half of trials 
  trllength = cat(1,data.time{:});
  trllength = trllength(:,2) - trllength(:,1);
  trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
  trlspikes = full(cat(2,trlspikes{:}));
  trlspikerate = bsxfun(@rdivide,trlspikes,trllength');
  selind = sum(trlspikerate>1,2)>round(numel(data.trial)/2);
  for itrial = 1:numel(data.trial)
    data.trial{itrial} = data.trial{itrial}(selind,:);
  end
  data.label = data.label(selind);
  

  
  % set fns
  fourfn = [info.savepath currsess '_' 'fourier' '_' fnnamesuffix '.mat'];
  datafn = [info.savepath currsess '_' 'dataetc' '_' fnnamesuffix '.mat'];
  if ~exist(fourfn,'file') || ~exist(datafn,'file')
    %fourier = roe_fourierspiketrain_direct(data,freqoi,timwin,true);
    fourier = roe_fourierspiketrain_conv(data,freqoi,timwin,'segmentsparseconv');
    
    % save fourier and data
    save(fourfn,'fourier','-v7.3')
    
    % save data and stuff
    data = rmfield(data,'trial');
    save(datafn,'data','freqoi','timwin')
  end
   
   
  % get pruned data sets if desired
  prunemode = {'1of4','2of4','3of4','4of4'}; % {'odd','even'}  {'1of4','2of4','3of4','4of4'}
  for iprune = 1:numel(prunemode)
    % set fns
    currfnsuff = ['prunemode' '_' prunemode{iprune}];
    fourprfn = [info.savepath currsess '_' 'fourier' '_' fnnamesuffix '_' currfnsuff '.mat'];
    dataprfn = [info.savepath currsess '_' 'dataetc' '_' fnnamesuffix '_' currfnsuff '.mat'];
    if ~exist(fourprfn,'file') || ~exist(dataprfn,'file')
      
      % read in data
      data = rmr_crcnsssc3_readspike([info.datapath currsess],triallength);
      
      % select units with firing rate at least 1Hz in half of trials
      trllength = cat(1,data.time{:});
      trllength = trllength(:,2) - trllength(:,1);
      trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
      trlspikes = full(cat(2,trlspikes{:}));
      trlspikerate = bsxfun(@rdivide,trlspikes,trllength');
      selind = sum(trlspikerate>1,2)>round(numel(data.trial)/2);
      for itrial = 1:numel(data.trial)
        data.trial{itrial} = data.trial{itrial}(selind,:);
      end
      data.label = data.label(selind);
      
      % prune
      datapr = rmr_prunespiketrain(data,prunemode{iprune});
      
      % get fourier
      %fourier = roe_fourierspiketrain_direct(datapr,freqoi,timwin,true);
      fourier = roe_fourierspiketrain_conv(datapr,freqoi,timwin,'segmentsparseconv');

      
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





