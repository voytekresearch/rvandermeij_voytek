function rmr_crcnsssc3_getcorrgram



% get info
info = rmr_crcnsssc3_info;
info.datapath = '/Users/roemervandermeij/Work/Data/tmpdata/CRCNS_SSC3/';

%
triallength = 50;
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
  trlspikerate = trlspikes ./ trllength.';
  selind = sum(trlspikerate>1,2)>round(numel(data.trial)/2);
  for itrial = 1:numel(data.trial)
    data.trial{itrial} = data.trial{itrial}(selind,:);
  end
  data.label = data.label(selind);

    
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
    fn = [info.savepath currsess '_' 'contcorrgram' '_' timename '_' 'keeptrials' '.mat'];
  else
    fn = [info.savepath currsess '_' 'contcorrgram' '_' timename '.mat'];
  end
  save(fn,'contcorrgram','time','data','-v7.3');
end









function playground






% get info
info = rmr_crcnsssc3_info;
info.datapath = '/Users/roemervandermeij/Work/Data/tmpdata/CRCNS_SSC3/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/crcnsssc3/';

%
triallength = 50;
for    isess = 11%     :numel(info.session) % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  
  % read in data
  data = rmr_crcnsssc3_readspike([info.datapath currsess],triallength);
  
  % select units with firing rate at least 1Hz in half of trials 
  trllength = cat(1,data.time{:});
  trllength = trllength(:,2) - trllength(:,1);
  trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
  trlspikes = full(cat(2,trlspikes{:}));
  trlspikerate = trlspikes ./ trllength.';
  selind = sum(trlspikerate>1,2)>round(numel(data.trial)/2);
  for itrial = 1:numel(data.trial)
    data.trial{itrial} = data.trial{itrial}(selind,:);
  end
  data.label = data.label(selind);
  
  % load
  %fn = [info.savepath currsess '_' 'contcorrgram' '_' '-20to20ms' '.mat'];
  fn = [info.savepath currsess '_' 'contcorrgram' '_' '-20to20ms_keeptrials' '.mat'];
  load(fn)
  
  % average over trials
  ntrial = size(contcorrgram,1);
  contcorrgramall = squeeze(sum(contcorrgram,1));
  contcorrgramc1  = squeeze(sum(contcorrgram(1:(ntrial/2),:,:,:),1));
  contcorrgramc2  = squeeze(sum(contcorrgram((ntrial/2)+1:end,:,:,:),1));
  
  % set
  nunit = size(contcorrgramall,1);
  nbins = size(contcorrgramall,3);
  
  % eliminate auto corr grams
  for iunit = 1:nunit
    contcorrgramall(iunit,iunit,:) = NaN(1,nbins);
    contcorrgramc1(iunit,iunit,:)  = NaN(1,nbins);
    contcorrgramc2(iunit,iunit,:)  = NaN(1,nbins);
  end
  
  % create labelcmb
  labelcmb = ft_channelcombination('all',data.label,1,2);
  
  % unfold
  contcorrgramall = permute(contcorrgramall,[3 1 2]); % bins, seedunit targunit
  contcorrgramall = reshape(contcorrgramall,[nbins nunit.^2]); % unfolded over channel pairs
  contcorrgramc1  = permute(contcorrgramc1,[3 1 2]); % bins, seedunit targunit
  contcorrgramc1  = reshape(contcorrgramc1,[nbins nunit.^2]); % unfolded over channel pairs
  contcorrgramc2  = permute(contcorrgramc2,[3 1 2]); % bins, seedunit targunit
  contcorrgramc2  = reshape(contcorrgramc2,[nbins nunit.^2]); % unfolded over channel pairs
  
  
  % remove auto correlations
  remind = find(strcmp(labelcmb(:,1),labelcmb(:,2)));
  labelcmb(remind,:) = [];
  contcorrgramall(:,remind) = [];
  contcorrgramc1(:,remind) = [];
  contcorrgramc2(:,remind) = [];
  % remove based on cutoff of 25 min spikes
  remind = max(contcorrgramall,[],1)<25;
  labelcmb(remind,:) = [];
  contcorrgramall(:,remind) = [];
  contcorrgramc1(:,remind)  = [];
  contcorrgramc2(:,remind)  = [];
  % sort by maximum relative distance to median
  tmpcorrgram = sort(contcorrgramall,1);
  %maxreldev = (tmpcorrgram(end,:) ./ mean(tmpcorrgram,1)) ./ (mean(tmpcorrgram(end-7:end-2,:),1) ./ mean(tmpcorrgram,1));
  maxreldev = tmpcorrgram(end,:) ./ mean(tmpcorrgram(end-7:end-2,:),1);
  [maxreldev sortind] = sort(maxreldev,'descend');
  contcorrgramall = contcorrgramall(:,sortind);
  contcorrgramc1  = contcorrgramc1(:,sortind);
  contcorrgramc2  = contcorrgramc2(:,sortind);
  labelcmb = labelcmb(sortind,:);
  
  
  % plot the highest channel pairs
  nsel = 50;
  for iset = 1:5%ceil(size(contcorrgramall,2)/nsel)
    figure('numbertitle','off','name',[currsess ' nz-pairs: ' num2str(sum(sum(contcorrgramall,1)~=0)) ' set: ' num2str(iset) ' nunit: ' num2str(numel(data.label))]);
    for ipair = 1:nsel
      currind = (1:nsel) + ((iset-1)*nsel);
      subplot(ceil(sqrt(nsel)),ceil(sqrt(nsel)),ipair)
      plot(time*1000,contcorrgramc1(:,currind(ipair)))
      hold on
      plot(time*1000,contcorrgramc2(:,currind(ipair)),'r')
      %set(gca,'ylim',[0 maxcount])
      title([labelcmb{currind(ipair),1} '   ' labelcmb{currind(ipair),2}],'interpreter','none')
    end
  end
  
end






