function rmr_crcnspvc8_getfourier


    
    
% get info
info = rmr_crcnspvc8_info;
%info.datapath = '/Users/roemer/Work/Data/tmpdata/CRCNS_PVC8/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/crcnspvc8/';

% set name suffix
fnnamesuffix = 'four_conv_timeavg_timwin20ms_freq50-50-900_first270natonly';
%fnnamesuffix = 'four_conv_timeavg_timwin50ms_freq20-20-400_first270natonly';
fnnamesuffix = 'four_conv_timeavg_timwin20ms_freq50-50-900_first20natonly';

% fake fourier settings
freqoi  = 50:50:900;
timwin  = .020;
%freqoi  = 20:20:400;
%timwin  = .050;

% read data and get fourier
% nprune    = 20;
% pruneperc = 10;
trialdeftype = 'stimonly';
for    isess = 1%      :numel(info.session); % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  
  % read in data
  fsample = 20000; 
  data = rmr_crcnspvc8_readspike([info.datapath currsess],fsample,trialdeftype);
  
  % merge repetitions
  datanew = data;
  datanew.trial = cell(1,max(data.trialinfo(:,1)));
  datanew.time  = cell(1,max(data.trialinfo(:,1)));
  datanew.sampleinfo = [];
  datanew.trialinfo = data.trialinfo(data.trialinfo(:,2)==1,:);
  padding = sparse(zeros(numel(data.label),0.1*data.fsample));
  for iimg = 1:max(data.trialinfo(:,1))
    ind = data.trialinfo(:,1)==iimg;
    currtrial = data.trial(ind);
    for irep = 1:numel(currtrial)
      currtrial{irep} = [currtrial{irep} padding];
    end
    datanew.trial{iimg} = cat(2,currtrial{:});
    datanew.time{iimg}  = [0 size(datanew.trial{iimg},2)./data.fsample];
  end
  data = datanew;
  
  
  nsel = 20;
  trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
  trlspikes = sum(full(cat(2,trlspikes{:})),1);
  imgspikes = sum(reshape(trlspikes,[20 numel(trlspikes)./20]),1);
  [dum,ind] = sort(imgspikes,'descend');
  imgsel    = [ind(find(data.trialinfo(1+((ind-1)*20),4)==1,nsel/2)) ind(find(data.trialinfo(1+((ind-1)*20),4)==2,nsel/2))];
  trialind  = bsxfun(@plus,[1:20]',((imgsel-1)*20));
  trialind  = sort(trialind(:));
  data.trial = data.trial(trialind);
  data.time  = data.time(trialind);
  data.sampleinfo  = data.sampleinfo(trialind,:);
  data.trialinfo   = data.trialinfo(trialind,:);

  % set fns
  fourfn = [info.savepath currsess '_' trialdeftype '_' 'fourier' '_' fnnamesuffix '.mat'];
  datafn = [info.savepath currsess '_' trialdeftype '_' 'dataetc' '_' fnnamesuffix '.mat'];
  if ~exist(fourfn,'file') || ~exist(datafn,'file')
    %fourier = roe_fourierspiketrain_direct(data,freqoi,timwin,true);
    fourier = roe_fourierspiketrain_conv(data,freqoi,timwin,'sparseconv');
    
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
    fourprfn = [info.savepath currsess '_' trialdeftype '_' 'fourier' '_' fnnamesuffix '_' currfnsuff '.mat'];
    dataprfn = [info.savepath currsess '_' trialdeftype '_' 'dataetc' '_' fnnamesuffix '_' currfnsuff '.mat'];
    if ~exist(fourprfn,'file') || ~exist(dataprfn,'file')
      
      % fetch data again
      fsample = 20000;
      data = rmr_crcnspvc8_readspike([info.datapath currsess],fsample,trialdeftype);
      
      % select based on most active small/big images
      trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
      trlspikes = sum(full(cat(2,trlspikes{:})),1);
      imgspikes = sum(reshape(trlspikes,[20 numel(trlspikes)./20]),1);
      [dum,ind] = sort(imgspikes,'descend');
      imgsel    = [ind(find(data.trialinfo(1+((ind-1)*20),4)==1,nsel/2)) ind(find(data.trialinfo(1+((ind-1)*20),4)==2,nsel/2))];
      trialind  = bsxfun(@plus,[1:20]',((imgsel-1)*20));
      trialind  = sort(trialind(:));
      data.trial = data.trial(trialind);
      data.time  = data.time(trialind);
      data.sampleinfo  = data.sampleinfo(trialind,:);
      data.trialinfo   = data.trialinfo(trialind,:);
      
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
end












%   % select based on most active small/big images
%   nsel = 20;
%   trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
%   trlspikes = sum(full(cat(2,trlspikes{:})),1);
%   imgspikes = sum(reshape(trlspikes,[20 numel(trlspikes)./20]),1);
%   [dum,ind] = sort(imgspikes,'descend');
%   imgsel    = [ind(find(data.trialinfo(1+((ind-1)*20),4)==1,nsel/2)) ind(find(data.trialinfo(1+((ind-1)*20),4)==2,nsel/2))];
%   trialind  = bsxfun(@plus,[1:20]',((imgsel-1)*20));
%   trialind  = sort(trialind(:));
%   data.trial = data.trial(trialind);
%   data.time  = data.time(trialind);
%   data.sampleinfo  = data.sampleinfo(trialind,:);
%   data.trialinfo   = data.trialinfo(trialind,:);

