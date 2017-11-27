function rmr_crcnspvc8_getcorrgram



    
% get info
info = rmr_crcnspvc8_info;
%info.datapath = '/Users/roemer/Work/Data/tmpdata/CRCNS_PVC8/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/crcnspvc8/';

%
fnnamesuffix = 'first20natonly';

% 
trialdeftype = 'stimonly';
for    isess = 1  %    :numel(info.session); % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  
  % read in data
  fsample = 20000; 
  data = rmr_crcnspvc8_readspike([info.datapath currsess],fsample,trialdeftype);
    
  % select based on most active small/big images
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
  
  % get corrgram
  cfg = [];
  cfg.range      = [-20 20] ./ 1000; % in s
  cfg.keeptrials = 'yes';
  cfg.gaussfwhm  = 0.5 ./ 1000;
  cfg.fsample    = 4000;
  [contcorrgram,time] = rmr_contcorrgramspiketrain(cfg,data);
  
  % save
  data = rmfield(data,'trial');
  timename = [num2str(time(1)*1000) 'to' num2str(time(end)*1000) 'ms'];
  if strcmp(cfg.keeptrials,'yes')
    fn = [info.savepath currsess '_' trialdeftype '_' 'contcorrgram' '_' timename '_' fnnamesuffix '_' 'keeptrials' '.mat'];
  else
    fn = [info.savepath currsess '_' trialdeftype '_' 'contcorrgram' '_' timename '_' fnnamesuffix '.mat'];
  end
  save(fn,'contcorrgram','time','data','-v7.3');
end








