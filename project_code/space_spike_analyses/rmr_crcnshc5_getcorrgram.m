function rmr_crcnshc5_getcorrgram



% get info
info = rmr_crcnshc5_info;
info.datapath = '/Users/roemervandermeij/Work/Data/tmpdata/CRCNS_HC5/';

% 
segsel = 1;
trialdeftype = 'mazerun';%'5secseg';%'mazerun'
for    isess = 1  %    :numel(info.session); % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  
  % read in data
  fsample = 20000; % oversample intentionally
  data = rmr_crcnshc5_readspike([info.datapath currsess],fsample,'mazerun',0,1,segsel);
  
  % select units with firing rate above 0.2 hz
  trllength = cat(1,data.time{:});
  trllength = trllength(:,2) - trllength(:,1);
  trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
  trlspikes = full(cat(2,trlspikes{:}));
  spikerate = sum(trlspikes,2) ./ sum(trllength);
  selind = spikerate>0.2;
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
    fn = [info.savepath currsess '_' trialdeftype '_' 'contcorrgram' '_' timename '_' num2str(segsel) 'segsel' '_' 'keeptrials' '.mat'];
  else
    fn = [info.savepath currsess '_' trialdeftype '_' 'contcorrgram' '_' timename '_' num2str(segsel) 'segsel' '.mat'];
  end
  save(fn,'contcorrgram','time','data','-v7.3');
end









function playground






% get info
info = rmr_crcnshc5_info;
info.datapath = '/Users/roemer/Work/Data/tmpdata/CRCNS_HC5/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/crcnshc5/';

% 
trialdeftype = 'mazerun';%'5secseg';%'mazerun'
for    isess = 1     :numel(info.session); % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  
  %   % read in data
  fsample = 20000; % oversample intentionally
  data = rmr_crcnshc5_readspike([info.datapath currsess],fsample,'mazerun',0,1,segsel);

  % select units with firing rate above 0.2 hz
  trllength = cat(1,data.time{:});
  trllength = trllength(:,2) - trllength(:,1);
  trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
  trlspikes = full(cat(2,trlspikes{:}));
  spikerate = sum(trlspikes,2) ./ sum(trllength);
  selind = spikerate>0.2;
  for itrial = 1:numel(data.trial)
    data.trial{itrial} = data.trial{itrial}(selind,:);
  end
  data.label = data.label(selind);
  
  
  % load
  fn = [info.savepath currsess '_' 'corrgram' '_' trialdeftype '_' '-20msto20ms1msbins_3segsel_keeptrials' '.mat'];
  %fn = [info.savepath currsess '_' 'corrgram' '_' trialdeftype '_' '-25msto25ms0.1msbins_3segsel_keeptrials' '.mat'];
  %fn = [info.savepath currsess '_' 'corrgram' '_' trialdeftype '_' '-100msto100ms1msbins_3segsel_keeptrials' '.mat'];
  load(fn)
  
  % average over trials
  corrgram = squeeze(sum(corrgram,1));
  
  % set
  nunit = size(corrgram,1);
  nbins = size(corrgram,3);
  
  % eliminate auto corr grams
  for iunit = 1:nunit
    corrgram(iunit,iunit,:) = NaN(1,nbins);
  end
  
  % create labelcmb
  labelcmb = ft_channelcombination('all',data.label,1,2);
  
  % unfold
  corrgram = permute(corrgram,[3 1 2]); % bins, seedunit targunit
  corrgram = reshape(corrgram,[nbins nunit.^2]); % unfolded over channel pairs 
  
  % remove auto correlations
  remind = find(strcmp(labelcmb(:,1),labelcmb(:,2)));
  labelcmb(remind,:) = [];
  corrgram(:,remind) = [];
  % remove based on cutoff of 25 min spikes
  remind = max(corrgram,[],1)<25;
  labelcmb(remind,:) = [];
  corrgram(:,remind) = [];
  % sort by maximum relative distance to median
  tmpcorrgram = sort(corrgram,1);
  %maxreldev = (tmpcorrgram(end,:) ./ mean(tmpcorrgram,1)) ./ (mean(tmpcorrgram(end-7:end-2,:),1) ./ mean(tmpcorrgram,1));
  maxreldev = tmpcorrgram(end,:) ./ mean(tmpcorrgram(end-7:end-2,:),1);
  [maxreldev sortind] = sort(maxreldev,'descend');
  corrgram = corrgram(:,sortind);
  labelcmb = labelcmb(sortind,:);

  
  % plot the highest channel pairs
  nsel = 50;
  maxcount = max(corrgram(:));
  for iset = 1%:ceil(size(corrgram,2)/nsel)
    figure('numbertitle','off','name',[currsess ' nz-pairs: ' num2str(sum(sum(corrgram,1)~=0)) ' set: ' num2str(iset)]);
    for ipair = 1:nsel
      currind = (1:nsel) + ((iset-1)*nsel);
      subplot(ceil(sqrt(nsel)),ceil(sqrt(nsel)),ipair)
      plot(time*1000,corrgram(:,currind(ipair)))
      %set(gca,'ylim',[0 maxcount])
      title([labelcmb{currind(ipair),1} '   ' labelcmb{currind(ipair),2}],'interpreter','none')
    end
  end
  
end





% GG.069
% unit50_11 & unit53_11
seedind = find(strcmp(data.label,'unit50_11'));
targind = find(strcmp(data.label,'unit53_11'));
seedind = 16;
targind = 1;
ntrial  = numel(data.time);
figure('numbertitle','off','name',[currsess ' ' data.label{seedind} ' ' data.label{targind}]);
maxcount = round(max(max(squeeze(corrgram(:,seedind,targind,:))))*1.2);
for itrial = 1:ntrial
  currcg = squeeze(corrgram(itrial,seedind,targind,:));
  subplot(ceil(sqrt(ntrial)),ceil(sqrt(ntrial)),itrial)
  bar(timebins,currcg,1)
  set(gca,'ylim',[0 maxcount])
  set(gca,'xlim',[timebins(1) timebins(end)])
  title(['trial' num2str(itrial)])
end
subplot(ceil(sqrt(ntrial)),ceil(sqrt(ntrial)),itrial+1)
currcg = squeeze(sum(corrgram(:,seedind,targind,:),1));

bar(timebins*1000,currcg,1)
set(gca,'xtick',binedges(1:round(numel(timebins)/10):end)*1000)
xlabel(['time distance from spikes of unit' num2str(seedind) ' (ms)'])
ylabel(['counts of unit' num2str(targind)])
title(['correlogram: seed = unit' num2str(seedind) ', target = unit' num2str(targind)])




figure
currcg = squeeze(sum(corrgram(:,seedind,targind,:),1));
bar(timebins,currcg,1)
set(gca,'xlim',[timebins(1) timebins(end)])
title('all trials')




figure('numbertitle','off','name',[currsess ' ' data.label{seedind} ' ' data.label{targind}]);
axh   = [];
axlim = [];
for itrial = 1:size(fourier,3)
  subplot(ceil(sqrt(ntrial)),ceil(sqrt(ntrial)),itrial)
  csd = zeros(size(fourier,1),size(fourier,1),size(fourier,2));
  for ifreq = 1:numel(freqoi)
    tmpfour = squeeze(fourier(:,ifreq,itrial,:));
    tmpfour(:,isnan(tmpfour(1,:))) = [];
    csd(:,:,ifreq) = csd(:,:,ifreq) + ((tmpfour * tmpfour')./size(fourier,3));
  end
  curr = squeeze(csd(seedind,targind,:));
  roe_compass(curr);
  title(['trial' num2str(itrial)])
  axh(itrial)     = gca;
  axlim(itrial,:) = axis;
end
% equalize axis
[dum maxind] = max(axlim(:,4));
for itrial = 1:size(fourier,3)
  axis(axh(itrial),axlim(maxind,:));
end
csd = zeros(size(fourier,1),size(fourier,1),size(fourier,2));
for itrial = 1:size(fourier,3)
  for ifreq = 1:numel(freqoi)
    tmpfour = squeeze(fourier(:,ifreq,itrial,:));
    tmpfour(:,isnan(tmpfour(1,:))) = [];
    csd(:,:,ifreq) = csd(:,:,ifreq) + ((tmpfour * tmpfour')./size(fourier,3));
  end
end
subplot(ceil(sqrt(ntrial)),ceil(sqrt(ntrial)),itrial+1)
curr = squeeze(csd(seedind,targind,:));
roe_compass(curr);
title('all trials')
  
%     'unit39_10'    'unit33_9' 
%     'unit34_9'     'unit22_9' 
%     'unit53_11'    'unit50_11'
%     'unit38_10'    'unit34_9' 
%     'unit50_11'    'unit49_11'

ind = match_str(data.label,{'unit33_9'    'unit50_11'    'unit58_12'    'unit28_9'    'unit9_6'    'unit23_9'    'unit51_11'    'unit15_7'    'unit43_11'    'unit17_8'});












for icomp = comporder
  figure('numbertitle','off','name',['comp' num2str(icomp)]);
  
  [dum sortind] = sort(A(:,icomp),'descend');
  chanind = sortind(1:4);
  
  subplot(3,3,1)
  plot(1:nchan,A(:,icomp));
  ylim = get(gca,'ylim');
  ylim = round([0 ylim(2)*1.05]*10)./10;
  set(gca,'ylim',ylim,'xlim',[1 nchan],'ytick',ylim)
  xlabel('units')
  ylabel('loading (au)')
  title('spatial map of network')
  subplot(3,3,2)
  % C t-test: enter target or not
  subplot(3,3,2)
  avgdist = cat(2,mean(nwaycomp.trialinfo(:,end-1:end),2),mean(nwaycomp.trialinfo(:,end-1:end),2));
  avgdist(nwaycomp.trialinfo(:,3)==1,1) = NaN;
  avgdist(nwaycomp.trialinfo(:,3)==2,2) = NaN;
  [ax h1 h2] = plotyy(1:ntrial,C(:,icomp),1:ntrial,avgdist);
  ylim = get(ax(1),'ylim');
  ylim = [0 ylim(2)*1.05];
  set(ax(1),'ylim',ylim,'xlim',[1 ntrial],'ytick',ylim)
  ylim = get(ax(2),'ylim');
  ylim = [0 ylim(2)*1.05];
  set(ax(2),'ylim',ylim,'xlim',[1 ntrial],'ytick',ylim)
  set(h2(1),'marker','o')
  set(h2(2),'marker','o')
  set(h1,'marker','d')
  xlabel('trial')
  ylabel('loading (au)')
  ind1 = nwaycomp.trialinfo(:,3)==1;
  ind2 = ~ind1;
  m1 = mean(C(ind1,icomp));
  m2 = mean(C(ind2,icomp));
  sem1 = std(C(ind1,icomp)) ./ sqrt(sum(ind1));
  sem2 = std(C(ind2,icomp)) ./ sqrt(sum(ind2));
  [h p,ci,stats] = ttest2(C(ind1,icomp),C(ind2,icomp));
  title(['network strength right/left trials, mean+sem, t = ' num2str(stats.tstat,'%1.2f') ',  corr pos ' num2str(corr(C(:,icomp),mean(nwaycomp.trialinfo(:,end-1:end),2)),'%1.2f')])
  % plot S
  subplot(3,3,3)
  hold on
  [dum sortind] = sort(A(:,icomp),'descend');
  [ax h1 h2] = plotyy(1:5,(S(sortind(1:5),icomp)-S(sortind(1),icomp))*1000,1.5:4.5,diff((S(sortind(1:5),icomp)-S(sortind(1),icomp))*1000));
  set(ax(1),'ylim',[-.010 .010]*1000,'xlim',[1 5],'xtick',1:5)
  set(ax(1),'ytick',(-.01 :0.005: .01)*1000)
  set(ax(2),'ylim',[-.005 .005]*1000,'xlim',[1 5],'xtick',1:5)
  set(ax(2),'ytick',(-.005 :0.0025: .005)*1000)
  set(h2,'linestyle','--','marker','o','markerfacecolor', [0.85 0.325 0.098])
  xlabel('units sorted by spatial map strength')
  ylabel(ax(1),'time delay (ms)')
  ylabel(ax(2),'diff of time delay (ms)')
  title('time delays of units')
  
  count = 0;
  for iseed = 1:numel(chanind)
    for itarg = (iseed+1):numel(chanind)
      count = count+1;
      subplot(3,3,3+count);
      currcg = squeeze(sum(corrgram(:,chanind(iseed),chanind(itarg),:),1));
      h = bar(timebins*1000,currcg,1);
      set(h,'facecolor',[0 .75 1])
      set(gca,'xtick',binedges(1:round(numel(timebins)/4):end)*1000)
      ylim = get(gca,'ylim');
      ylim = round([0 ylim(2)*1.05]);
      set(gca,'ylim',ylim,'ytick',ylim)
      xlabel(['time from seed unit spikes (ms)'])
      ylabel(['target unit spike counts'])
      title(['correlogram: seed = unit' num2str(iseed) ', target = unit' num2str(itarg)])
    end
  end
end









for icomp = comporder
  figure('numbertitle','off','name',['comp' num2str(icomp)]);
  [dum sortind] = sort(A(:,icomp),'descend');
  chanind = sortind(1:4);
  
  count = 0;
  ind1 = nwaycomp.trialinfo(:,3)==1;
  ind2 = ~ind1;
  for iseed = 1:numel(chanind)
    for itarg = (iseed+1):numel(chanind)
      h = [];
      for iset = 1:2
        if iset==1
          currind = ind1;
        else
          currind = ind2;
        end
        count = count+1;
        subplot(4,3,count);
        currcg = squeeze(corrgram(currind,chanind(iseed),chanind(itarg),:));
        trialtime = cellfun(@max,data.time(currind));
        currcg = currcg ./ repmat(trialtime',[1 size(currcg,2)]);
        currcg = mean(currcg,1);
        %currcg = sum(currcg,1);
        hb = bar(timebins*1000,currcg,1);
        h(iset) = gca;
        set(hb,'facecolor',[0 .75 1])
        set(gca,'xtick',binedges(1:round(numel(timebins)/4):end)*1000)
        ylim = get(gca,'ylim');
        ylim = round([0 ylim(2)*1.05]*100)./100;
        set(gca,'ylim',ylim,'ytick',ylim)
        xlabel(['time from seed unit spikes (ms)'])
        ylabel(['avg target unit spikes/s'])
        title(['correlogram: seed = unit' num2str(iseed) ', target = unit' num2str(itarg)])
      end
      ylim = [0 round(max(cellfun(@max,get(h,'ylim')))*100)./100];
      set(h,'ylim',ylim,'ytick',ylim);
    end
  end
end






















