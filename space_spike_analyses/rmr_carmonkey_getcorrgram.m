function rmr_carmonkey_getcorrgram



% get info
info = rmr_carmonkey_info;
info.datapath = '/Users/roemervandermeij/Work/Data/tmpdata/carmena_monkey/';

% set name suffix
fnnamesuffix = '1hzfiring_gotilltar';
%fnnamesuffix = '1hzfiring_centertillgo';


% extract trialdeftype
ind = strfind(fnnamesuffix,'_');
trialdeftype = fnnamesuffix(ind(end)+1:end);

% read data and get fourier, 1 session
keeptrials = 'yes';
for      isubj = 1     :numel(info.subj)
  for    isess = 1     :numel(info.session.(info.subj{isubj}))
    % set currs
    currsubj = info.subj{isubj};
    currsess = info.session.(currsubj){isess};
    
    % set corrgram options
    cfg = [];
    cfg.range      = [-20 20] ./ 1000; % in s
    cfg.keeptrials = 'yes';
    cfg.gaussfwhm  = 0.5 ./ 1000;
    cfg.fsample    = 4000;
    timename = [num2str(cfg.range(1)*1000) 'to' num2str(cfg.range(end)*1000) 'ms'];
    
    % check whether done
    if strcmp(cfg.keeptrials,'yes')
      fn = [info.savepath currsess '_' fnnamesuffix '_' 'contcorrgram' '_' timename '_' 'keeptrials' '.mat'];
    else
      fn = [info.savepath currsess '_' fnnamesuffix '_' 'contcorrgram' '_' timename  '.mat'];
    end
    if exist(fn,'file')
      continue
    end
    

    % fetch data
    fsample = 40000; % oversample intentionally
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
    
    
    % get corrgram
    [contcorrgram,time] = rmr_contcorrgramspiketrain(cfg,data);
    
    % save
    data = rmfield(data,'trial');
    save(fn,'contcorrgram','time','data','-v7.3');
  end
end
















function playground




% get info
info = rmr_carmonkey_info;
%info.datapath = '/Users/roemer/Work/Data/tmpdata/carmena_monkey/';

% set name suffix
fnnamesuffix = '1hzfiring_gotilltar';
%fnnamesuffix = '1hzfiring_centertillgo';

% extract trialdeftype
ind = strfind(fnnamesuffix,'_');
trialdeftype = fnnamesuffix(ind(end)+1:end);

% read data and get fourier, 1 session
for      isubj = 1     :numel(info.subj)
  for    isess = 1     ; % 1 session for now
    
    % set currs
    currsubj = info.subj{isubj};
    currsess = info.session.(currsubj){isess};
    
    
    % load
    %fn = [info.savepath currsess '_' 'corrgram' '_' fnnamesuffix '_' '-100msto100ms1msbins_keeptrials' '.mat'];
    fn = [info.savepath currsess '_' 'corrgram' '_' fnnamesuffix '_' '-20msto20ms1msbins_keeptrials' '.mat'];
    %fn = [info.savepath currsess '_' 'corrgram' '_' fnnamesuffix '_' '-20msto20ms1msbins' '.mat'];
    %fn = [info.savepath currsess '_' 'corrgram' '_' fnnamesuffix '_' '-10msto10ms0.1msbins' '.mat'];
    load(fn)
    
    % average over trials
    corrgram = squeeze(sum(corrgram,1));
    
    % set
    nunit = size(corrgram,1);
    nbins = size(corrgram,3);
    
    % create labelcmb
    labelcmb = ft_channelcombination('all',data.label,1,2);
    
    % unfold
    corrgram = permute(corrgram,[3 1 2]); % bins, seedunit targunit
    corrgram = reshape(corrgram,[size(corrgram,1) size(corrgram,2) * size(corrgram,3)]); % unfolded over channel pairs
    
    % sort by something
    %[tmpsort, sortind] = sort(max(corrgram,[],1),'descend');
    %[tmpsort, sortind] = sort(sqrt(nansum(abs(corrgram ./ repmat(nansum(corrgram,1),[nbins 1 1])).^2,1)),'descend'); % sort by L2 norm after norming to sum 1
    %     [tmpsort, sortind] = sort(max(corrgram,[],1) ./ nanmean(corrgram,1),'descend');
    %     sortind = sortind([find(~isnan(tmpsort)) find(isnan(tmpsort))]);
    %     tmpsort = tmpsort([find(~isnan(tmpsort)) find(isnan(tmpsort))]);
    %     sortind = sortind([find(~isinf(tmpsort)) find(isinf(tmpsort))]);
    %     tmpsort = tmpsort([find(~isinf(tmpsort)) find(isinf(tmpsort))]);
    
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
    nsel = 120;
    for iset = 1%:5%ceil(size(corrgram,2)/nsel)
      figure('numbertitle','off','name',[currsess ' nz-pairs: ' num2str(sum(sum(corrgram,1)~=0)) ' set: ' num2str(iset)]);
      for ipair = 1:nsel
        currind = (1:nsel) + ((iset-1)*nsel);
        subplot(ceil(sqrt(nsel)),ceil(sqrt(nsel)),ipair)
        bar(timebins,corrgram(:,currind(ipair)),1)
        %set(gca,'ylim',[0 maxcount])
        title([labelcmb{currind(ipair),1} '   ' labelcmb{currind(ipair),2}],'interpreter','none')
      end
    end
    
    
  end
end













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
  ind1 = nwaycomp.trialinfo(:,1)~=12; % 7 == entering target, 12 == timeout
  ind2 = ~ind1;
  m1 = mean(C(ind1,icomp));
  m2 = mean(C(ind2,icomp));
  sem1 = std(C(ind1,icomp)) ./ sqrt(sum(ind1));
  sem2 = std(C(ind2,icomp)) ./ sqrt(sum(ind2));
  errorbar([m1 m2],[sem1 sem2])
  ylim = get(gca,'ylim');
  ylim = round([0 ylim(2)*1.1]*100)./100;
  set(gca,'ylim',ylim);
  set(gca,'xlim',[0.5 2.5])
  set(gca,'xtick',[1 2])
  set(gca,'ytick',ylim)
  set(gca,'xticklabel',{'cursor reached target (hit)','failed to reach target (miss)'})
  ylabel('loading (au)')
  [h p,ci,stats] = ttest2(C(ind1,icomp),C(ind2,icomp));
  title(['network strength hits/miss trials, mean+sem, t = ' num2str(stats.tstat,'%1.2f')])
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
  ind1 = nwaycomp.trialinfo(:,1)~=12;
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
        
        %         % plot as normal (mean) hist
        %         currcg = squeeze(corrgram(currind,chanind(iseed),chanind(itarg),:));
        %         currcg = mean(currcg,1);
        %         hb = bar(timebins*1000,currcg,1);
        
        %         % plot as spikes/sec
        %         currcg = squeeze(corrgram(currind,chanind(iseed),chanind(itarg),:));
        %         trialtime = cellfun(@max,data.time(currind));
        %         currcg = currcg ./ repmat(trialtime',[1 size(currcg,2)]);
        %         hb = bar(timebins*1000,nanmean(currcg,1),1);
        
        %         % plot as proportion of total trial-averaged spikes/sec in search window
        %         currcg = squeeze(corrgram(currind,chanind(iseed),chanind(itarg),:));
        %         trialtime = cellfun(@max,data.time(currind));
        %         currcg = currcg ./ repmat(trialtime',[1 size(currcg,2)]);
        %         currcg = nanmean(currcg,1);
        %         currcg = currcg ./ sum(currcg);
        %         hb = bar(timebins*1000,currcg,1);
        
        % plot as trial-averaged proportion of total  spikes in trial
        currcg = squeeze(corrgram(currind,chanind(iseed),chanind(itarg),:));
        currcg = bsxfun(@rdivide,currcg,sum(trlspikes([chanind(iseed),chanind(itarg)],currind),1)');
        currcg = nanmean(currcg,1);
        hb = bar(timebins*1000,currcg,1);
        
        %         % plot as spikes/sec fraction of trial-averaged normalized firing rate
        %         currcg = squeeze(corrgram(currind,chanind(iseed),chanind(itarg),:));
        %         trialtime = cellfun(@max,data.time(currind));
        %         currcg = currcg ./ repmat(trialtime',[1 size(currcg,2)]);
        %         currcg = bsxfun(@rdivide,currcg,mean(spikerate([chanind(iseed),chanind(itarg)],currind),1)');
        %         currcg = nanmean(currcg,1);
        %         hb = bar(timebins*1000,currcg,1);
        
        h(iset) = gca;
        set(hb,'facecolor',[0 .75 1])
        set(gca,'xtick',binedges(1:round(numel(timebins)/4):end)*1000)
        xlabel(['time from seed unit spikes (ms)'])
        %ylabel(['avg target unit spikes/s'])
        %ylabel(['avg prop.'])
        title(['cond' num2str(iset) ': unit' num2str(iseed) '-' 'unit' num2str(itarg)])
      end
      ylim = [0 round(max(cellfun(@max,get(h,'ylim')))*1000)./1000];
      set(h,'ylim',ylim,'ytick',ylim);
    end
  end
end










% fetch data
fsample = 40000; % oversample intentionally
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
trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
trlspikes = full(cat(2,trlspikes{:}));
spikerate = bsxfun(@rdivide,trlspikes,trialtime);










