function rmr_predfaceval_getpsdandtfr





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
  fn{1,1} = [info.savepath currsubj '_' 'tfrpow_4to60hz_3cyc'     '_'         'CTI'   '.mat'];
  fn{2,1} = [info.savepath currsubj '_' 'tfrpow_60to120hz_3cyc'   '_'         'CTI'   '.mat'];
  fn{3,1} = [info.savepath currsubj '_' 'tfrpow_50to200hz_100ms'  '_'         'CTI'   '.mat'];
  fn{1,2} = [info.savepath currsubj '_' 'tfrpow_4to60hz_3cyc'     '_' 'postface500'   '.mat'];
  fn{2,2} = [info.savepath currsubj '_' 'tfrpow_60to120hz_3cyc'   '_' 'postface500'   '.mat'];
  fn{3,2} = [info.savepath currsubj '_' 'tfrpow_50to200hz_100ms'  '_' 'postface500'   '.mat'];
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
        toilim   = [-0.3 1.2];
        minfreq  = 4;
        trialsel = {predtc==1, predtc==2};
      case {2} % postface
        toilim   = [1 1.9];
        minfreq  = 4;
        trialsel = {valencetc==1, valencetc==2, valencetc==1 & predtc==1, valencetc==1 & predtc==2, valencetc==2 & predtc==1, valencetc==2 & predtc==2};
    end
    
    % prep settings
    freq    = {minfreq:1:60     ,60:2:120           ,50:5:200                    };
    timwin  = {3./freq{1}       ,3./freq{2}         ,ones(size(freq{3})) .* .050};
    
    
    % all TFRs
    for iset = 1:3
      if ~exist(fn{iset,iperiod},'file')
        freqdata = [];
        for itrialset = 1:numel(trialsel)
          % get TFR
          cfg = [];
          cfg.channel    = 'all';
          cfg.trials     = find(trialsel{itrialset});
          cfg.keeptrials = 'no';
          cfg.pad        = ceil(max(cellfun(@numel,data.time)) ./ data.fsample);
          cfg.method     = 'mtmconvol';
          cfg.output     = 'pow';
          cfg.taper      = 'hanning';
          cfg.foi        = freq{iset};
          cfg.t_ftimwin  = timwin{iset};
          cfg.toi        = toilim(1):(1/250):toilim(2);
          freqdata{itrialset} = ft_freqanalysis(cfg,data);
        end
        
        % save
        save(fn{iset,iperiod},'freqdata','-v7.3');
        clear freqdata
      end
    end % iset
    
  end % iperiod
  
end













function playground





% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';


for     isubj = 2:3%  :numel(info.subj)
  
  
  % set
  currsubj   = info.subj{isubj};
  fn = [];
  fn{1,1} = [info.savepath currsubj '_' 'tfrpow_4to60hz_3cyc'     '_'         'CTI'   '.mat'];
  fn{2,1} = [info.savepath currsubj '_' 'tfrpow_60to120hz_3cyc'   '_'         'CTI'   '.mat'];
  fn{3,1} = [info.savepath currsubj '_' 'tfrpow_50to200hz_100ms'  '_'         'CTI'   '.mat'];
  fn{1,2} = [info.savepath currsubj '_' 'tfrpow_4to60hz_3cyc'     '_' 'postface500'   '.mat'];
  fn{2,2} = [info.savepath currsubj '_' 'tfrpow_60to120hz_3cyc'   '_' 'postface500'   '.mat'];
  fn{3,2} = [info.savepath currsubj '_' 'tfrpow_50to200hz_100ms'  '_' 'postface500'   '.mat'];
  
  % for each of the freqs, plot the tfr
  for iperiod = 2
    switch iperiod
      case {1} % CTI
        toilim  = [0.2 1.2];
        condname = {'pred','unpred'};
      case {2} % postface
        toilim  = [1.2 1.8];
        condname = {'fear','neat','fearpred','fearunpred','neutpred','neutunpred'};
    end
    for iset = [1]
      
      % load that shit up
      load(fn{iset,iperiod})
      
      for ifreqset = 1:2%:numel(freqdata)
        currfreqdata = freqdata{ifreqset};
        % baseline by hand
        bslfd = load(fn{iset,1});
        bslfd1 = bslfd.freqdata{1};
        bslfd2 = bslfd.freqdata{2};
        ind = nearest(bslfd.freqdata{1}.time,0); % all time axis are identical after hfa analysis
        bslcorr1 = nanmean(bslfd1.powspctrm(:,:,1:ind),3) * numel(bslfd1.cfg.trials);
        bslcorr2 = nanmean(bslfd2.powspctrm(:,:,1:ind),3) * numel(bslfd2.cfg.trials);
        bslcorr = nansum(cat(3,bslcorr1,bslcorr1),3) ./ (numel(bslfd1.cfg.trials)+numel(bslfd2.cfg.trials));
        currfreqdata.powspctrm = bsxfun(@rdivide,currfreqdata.powspctrm,bslcorr);
        %currfreqdata.powspctrm = log10(currfreqdata.powspctrm);
        % set some things
        nchan = numel(currfreqdata.label);
        
        % get layout
        cfg = [];
        cfg.layout = 'ordered';
        tmp = currfreqdata;
        tmp.label = tmp.label;
        lay = ft_prepare_layout(cfg,tmp);
        
        % reorder channels in lay
        ind = [];
        type = {'RAM','LAM','AM'};
        for itype = 1:numel(type)
          ind = [ind; find(strncmp(lay.label,type{itype},numel(type{itype})))];
        end
        ind = [ind; setdiff(1:nchan,ind)'];
        lay.label  = lay.label(ind);
        %lay.pos    = lay.pos(ind,:);
        lay.width  = lay.width(ind);
        lay.height = lay.height(ind);
        
        % plot
        if iperiod == 1
          figname = [currsubj '-' 'CTI'];
        else
          figname = [currsubj '-' 'postface'];
        end
        if iset == 1
          figname = [figname '-' '4to60'];
        elseif iset == 2
          figname = [figname '-' '60to120'];
        else
          figname = [figname '-' '50to200'];
        end
        figname = [figname '-' condname{ifreqset}];
        figure('numbertitle','off','name',figname)
        
        % get clim
        clim = [0 0];
        for ichan = 1:nchan
          clim(2) = max(clim(2),max(max(abs(squeeze(currfreqdata.powspctrm(ichan,:,:))))));
        end
        clim = round(clim);
        clim = [0 2];
        
        % plot
        for ichan = 1:nchan
          subplot(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan)
          cfg = [];
          cfg.channel     = ichan;
          cfg.colorbar    = 'no';
          cfg.zlim        = clim;
          cfg.interactive = 'no';
          ft_singleplotTFR(cfg,currfreqdata);
          axis off
          if iperiod == 1
            line([0 0],ylim,'color','k')
            line([.2 .2],ylim,'color','k')
          else
            line([1.2 1.2],ylim,'color','k')
          end
        end
        
        
        
        
        
      end
    end % iset
    
  end % iperiod
  
  
end % isubj
%%%%%%%%%%%%%%
















% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

for     isubj = 1  :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  
  % check whether already done
  fn = [];
  fn{1,1} = [info.savepath currsubj '_' 'tfrpow_4to60hz_3cyc'     '_'         'CTI'   '.mat'];
  fn{2,1} = [info.savepath currsubj '_' 'tfrpow_60to120hz_3cyc'   '_'         'CTI'   '.mat'];
  fn{3,1} = [info.savepath currsubj '_' 'tfrpow_50to200hz_100ms'  '_'         'CTI'   '.mat'];
  fn{1,2} = [info.savepath currsubj '_' 'tfrpow_4to60hz_3cyc'     '_' 'postface500'   '.mat'];
  fn{2,2} = [info.savepath currsubj '_' 'tfrpow_60to120hz_3cyc'   '_' 'postface500'   '.mat'];
  fn{3,2} = [info.savepath currsubj '_' 'tfrpow_50to200hz_100ms'  '_' 'postface500'   '.mat'];
  
  % skip if needed
  if ~exist(fn{3},'file')
    %continue
  end
  
  % for each of the freqs, plot the tfr
  for iset = 1:numel(fn)
    
    % load that shit up
    load(fn{iset})
    
    % set some things
    freq  = freqdata.freq;
    nfreq = numel(freq);
    nchan = numel(freqdata.label);
    
    % change channel order to plot some before others
    ind = [];
    type = {'RAM','LAM','ROF','LOF'};
    for itype = 1:numel(type)
      ind = [ind; find(strncmp(freqdata.label,type{itype},3))];
    end
    ind = [ind; setdiff(1:nchan,ind)'];
    freqdata.label = freqdata.label(ind);
    freqdata.powspctrm = freqdata.powspctrm(:,ind,:,:);
    
    % fear/neut/pred/unpred/fearvsneut/predvsunpred
    for iplot = 1:6
      switch iplot
        case {1,2,3,4}
          switch iplot
            case 1
              trialind = freqdata.trialinfo(:,2) == 1; % fear
            case 2
              trialind = freqdata.trialinfo(:,2) == 2; % neut
            case 3
              trialind = freqdata.trialinfo(:,1) == 1; % pred
            case 4
              trialind = freqdata.trialinfo(:,1) == 2; % unpred
          end
          % mean
          cfg = [];
          cfg.trials = trialind;
          tmpfreq = ft_freqdescriptives(cfg,freqdata);
          % baseline
          cfg = [];
          cfg.baseline     = [-0.5 0];
          cfg.baselinetype = 'vssum';
          tmpfreq = ft_freqbaseline(cfg,tmpfreq);
        case 5 % fearvsneut
          % mean
          cfg = [];
          cfg.trials = freqdata.trialinfo(:,2) == 1; % fear
          tmpfreq1 = ft_freqdescriptives(cfg,freqdata);
          cfg.trials = freqdata.trialinfo(:,2) == 2; % neut
          tmpfreq2 = ft_freqdescriptives(cfg,freqdata);
          % norm
          tmpfreq = tmpfreq1;
          tmpfreq.powspctrm = (tmpfreq1.powspctrm-tmpfreq2.powspctrm) ./ (tmpfreq1.powspctrm+tmpfreq2.powspctrm);
        case 6 % predvsunpred
          % mean
          cfg = [];
          cfg.trials = freqdata.trialinfo(:,1) == 1; % pred
          tmpfreq1 = ft_freqdescriptives(cfg,freqdata);
          cfg.trials = freqdata.trialinfo(:,1) == 2; % unpred
          tmpfreq2 = ft_freqdescriptives(cfg,freqdata);
          % norm
          tmpfreq = tmpfreq1;
          tmpfreq.powspctrm = (tmpfreq1.powspctrm-tmpfreq2.powspctrm) ./ (tmpfreq1.powspctrm+tmpfreq2.powspctrm);
      end
      
      % plot
      figname = {'fear','neut','pred','unpred','fearvsneut','predvsunpred'};
      figure('numbertitle','off','name',[currsubj figname{iplot}])
      clim = [-max(abs(tmpfreq.powspctrm(:))) max(abs(tmpfreq.powspctrm(:)))] .* .5;
      for ichan = 1:nchan
        subplot(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan)
        cfg = [];
        cfg.channel     = ichan;
        cfg.colorbar    = 'no';
        cfg.zlim        = clim;
        cfg.interactive = 'no';
        if ismember(iplot,[1 2 5])
          cfg.xlim = [1 1.9];
        end
        ft_singleplotTFR(cfg,tmpfreq);
        set(gca,'ytick',[],'xtick',[0 1],'xticklabel',[],'ticklength',[0.3 0.3])
        if ichan == nchan
          set(gca,'ytick',[min(freq) min(freq)+((max(freq)-min(freq))/2) max(freq)])
          h = colorbar('eastoutside');
          set(h,'ticks',fix(clim*100)./100)
        end
      end
    end % iplot
    
    
    
    
  end % iset
end % isubj










































