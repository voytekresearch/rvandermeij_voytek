function rmr_predfaceval_gethfa





% fetch info
info = rmr_predfaceval_info;

%info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

for     isubj = 1  :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  
  % for artifacts been detected
  if ~info.(currsubj).trlartfctflg
    disp(['artifacts not yet detected for ' currsubj])
    continue
  end
  
  
  
  % check whether already done
  fn = [];
  fn{1} = [info.savepath currsubj '_' 'hfa100to200at50ms'     '_'         'CTI'   '.mat'];
  fn{2} = [info.savepath currsubj '_' 'hfa100to200at50ms'     '_'    'postface'   '.mat'];
  if all(cellfun(@exist,fn(:),repmat({'file'},[numel(fn) 1])))
    continue
  end
  
  % fetch and preprocess data
  cfg = [];
  cfg.demean      = 'yes';
  cfg.detrend     = 'yes';
  cfg.prestim     = 0.3;
  cfg.poststim    = 0.5+0.2;
  cfg.reref       = 'yes';
  cfg.refchannel  = 'all';
  cfg.bsfilter    = 'yes';
  cfg.lpfilter    = 'yes';
  data = rmr_predfaceval_fetchpreprocessdata(cfg,currsubj,info);
  
  
  % get trialdef
  valencetc = data.trialinfo(:,2); % fearful (1) or neutral (2) face trial
  predtc    = data.trialinfo(:,1); % pred (1) or unpred (2) trial
  
  % get it
  for iperiod = 1:2
    
    % combined over per condition and wplfs normalized within trials
    switch iperiod
      case {1} % CTI
        toilim   = [-0.2 1.2];
      case {2} % postface
        toilim   = [1 1.9];
    end
    
    % all TFRs
    if ~exist(fn{iperiod},'file')
      
      % get hfa
      freqoi = 100:10:200;
      timeoi = toilim(1):(1/250):toilim(2);
      timwin = 0.050;
      taper  = 'hanning';
      hfa = rmr_highfrequencyactivity(data,freqoi,timeoi,timwin,taper);
      
      % save
      save(fn{iperiod},'hfa','-v7.3');
      clear hfa
    end
    
  end % iperiod
  
  
end













function playground












% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

for     isubj = 1  :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  fn = [];
  fn{1} = [info.savepath currsubj '_' 'hfa100to200at50ms' '_' 'CTI' '.mat'];
  fn{2} = [info.savepath currsubj '_' 'hfa100to200at50ms' '_' 'postface' '.mat'];
  
  % skip if needed
  if ~exist(fn{2},'file')
    %continue
  end
  
  % for each of the freqs, plot the tfr
  for iperiod = 1:2
    % load that shit up
    load(fn{iperiod})
    
    % get trialdef
    valencetc = hfa.trialinfo(:,2); % fearful (1) or neutral (2) face trial
    predtc    = hfa.trialinfo(:,1); % pred (1) or unpred (2) trial
    
    switch iperiod
      case {1} % CTI
        trialsel = {predtc==1, predtc==2};
        setname  = {'predunpred'};
      case {2} % postface
        trialsel  = {valencetc==1, valencetc==2; valencetc==1 & predtc==1, valencetc==2 & predtc==1; valencetc==1 & predtc==2, valencetc==2 & predtc==2};
        setname  = {'fearneut','fearneutpred','fearneutunpred'};
    end
    % baseline by hand
    bslhfa = load(fn{1});
    bslhfa = bslhfa.hfa;
    ind = nearest(bslhfa.time{1},0); % all time axis are identical after hfa analysis
    bslcorr = cell(size(bslhfa.trial));
    for itrial = 1:numel(bslhfa.trial)
      bslcorr{itrial} = mean(bslhfa.trial{itrial}(:,1:ind),2);
    end
    for itrial = 1:numel(hfa.trial)
      hfa.trial{itrial} = hfa.trial{itrial} - repmat(bslcorr{itrial},[1 size(hfa.trial{itrial},2)]);
    end
    
    % do it per set
    for itrialset = 1:size(trialsel,1)
      % get timelocks
      timelock = [];
      cfg = [];
      cfg.trials = find(trialsel{itrialset,1});
      timelock{1} = ft_timelockanalysis(cfg,hfa);
      cfg.trials = find(trialsel{itrialset,2});
      timelock{2} = ft_timelockanalysis(cfg,hfa);
      
      % set some things
      nchan = numel(timelock{1}.label);
      
      % select channels
      ind = [];
      type = {'RAM','LAM','AM'};
      for itype = 1:numel(type)
        ind = [ind; find(strncmp(timelock{1}.label,type{itype},numel(type{itype})))];
      end
      ind = [ind; setdiff(1:nchan,ind)'];
      for iset = 1:numel(timelock)
        timelock{iset}.label = timelock{iset}.label(ind);
        timelock{iset}.avg   = timelock{iset}.avg(ind,:);
        timelock{iset}.var   = timelock{iset}.var(ind,:);
        timelock{iset}.dof   = timelock{iset}.dof(ind,:);
      end
      
      % plot in sets of two
      if iperiod == 1
        figure('numbertitle','off','name',[currsubj '-' setname{itrialset} '-' 'CTI'])
      elseif iperiod == 2
        figure('numbertitle','off','name',[currsubj '-' setname{itrialset} '-' 'postface'])
      end
      ylim = [min([timelock{1}.avg(:); timelock{2}.avg(:)]) max([timelock{1}.avg(:); timelock{2}.avg(:)])];
      ylim = round(ylim*10)./10;
      xlim = [timelock{1}.time(1) timelock{1}.time(end)];
      for ichan = 1:nchan
        subplot_tight(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan,0.02)
        hold on
        %
        nnanind = ~isnan(timelock{1}.avg(ichan,:));
        avg1 = timelock{1}.avg(ichan,nnanind);
        avg2 = timelock{2}.avg(ichan,nnanind);
        sem1 = sqrt(timelock{1}.var(ichan,nnanind)) ./ sqrt(mean(timelock{1}.dof(ichan,nnanind)));
        sem2 = sqrt(timelock{2}.var(ichan,nnanind)) ./ sqrt(mean(timelock{2}.dof(ichan,nnanind)));
        time = timelock{1}.time(nnanind);
        % plot patch, then lines
        x = [time time(end:-1:1)];
        y = [avg1-sem1 avg1(end:-1:1)+sem1(end:-1:1)];
        patch(x,y,'b','edgecolor','none','facealpha',0.15)
        y = [avg2-sem2 avg2(end:-1:1)+sem2(end:-1:1)];
        patch(x,y,'gr','edgecolor','none','facealpha',0.15)
        plot(time,avg1,'b')
        plot(time,avg2,'gr')
        title(timelock{1}.label{ichan})
        % plot axis
        %ylim = [-60 60];
        set(gca,'xlim',xlim,'ylim',ylim)
        if iperiod == 1
          line([0 0],ylim,'color','k')
          line([.2 .2],ylim,'color','k')
        else
          line([1.2 1.2],ylim,'color','k')
        end
        line(xlim,[0 0],'color','k')
        if ichan == nchan
          xlabel('time');ylabel('voltage')
          if iperiod == 1
            set(gca,'ytick',[ylim(1) 0 ylim(2)],'xtick',[xlim(1) 0 xlim(2)])
          else
            set(gca,'ytick',[ylim(1) 0 ylim(2)],'xtick',[xlim(1) 1.2 xlim(2)])
          end
        else
          set(gca,'ytick',[],'xtick',[])
        end
      end % nchan
      
      
    end % iset
    
  end % iperiod
  
  
end % isubj








































