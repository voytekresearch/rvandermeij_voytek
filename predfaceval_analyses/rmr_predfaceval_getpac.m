function rmr_predfaceval_getpac





% fetch info
info = rmr_predfaceval_info;

%info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

% design
centerampfreq = 'yes';
combtrials    = 'no';
for     isubj = 1:numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  
  % for artifacts been detected
  if ~info.(currsubj).trlartfctflg
    disp(['artifacts not yet detected for ' currsubj])
    continue
  end
  
  % check whether already done
  fn = [];
  fn{1,1} = [info.savepath currsubj '_' 'wplf_4to60'           '_'         'CTI' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials '.mat'];
  fn{2,1} = [info.savepath currsubj '_' 'wplf_4to60to70to130'  '_'         'CTI' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{3,1} = [info.savepath currsubj '_' 'wplf_4to60toHFA'      '_'         'CTI' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{1,2} = [info.savepath currsubj '_' 'wplf_6to60'           '_' 'postface500' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{2,2} = [info.savepath currsubj '_' 'wplf_6to60to70to130'  '_' 'postface500' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{3,2} = [info.savepath currsubj '_' 'wplf_6to60toHFA'      '_' 'postface500' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  if all(cellfun(@exist,fn(:),repmat({'file'},[numel(fn) 1])))
    continue
  end
  
  % combine data over sessions
  fndat = [info.savepath currsubj '_' 'allsess_bs540lp_detrend_0mspre_600mspost' '.mat'];
  if ~exist(fndat,'file')
    % fetch and preprocess data
    cfg = [];
    cfg.demean      = 'yes';
    cfg.detrend     = 'yes';
    cfg.prestim     = 0;
    cfg.poststim    = 0.5+0.1;
    cfg.reref       = 'yes';
    cfg.refchannel  = 'all';
    cfg.bsfilter    = 'yes';
    cfg.lpfilter    = 'yes';
    data = rmr_predfaceval_fetchpreprocessdata(cfg,currsubj,info);
    
    % save dat
    save(fndat,'data','-v7.3');
  else
    load(fndat)
  end
  
  
  % get trialdef
  valencetc = data.trialinfo(:,2); % fearful (1) or neutral (2) face trial
  predtc    = data.trialinfo(:,1); % pred (1) or unpred (2) trial
  
  % get PAC
  for iperiod = 1:2
    % combined over per condition and wplfs normalized within trials
    switch iperiod
      case {1} % CTI
        toilim  = [0.2 1.2];
        minfreq = 4;
        if istrue(combtrials)
          trialsel = ones(size(predtc));
        else
          trialsel = {predtc==1, predtc==2};
        end
      case {2} % postface
        toilim  = [1.2 1.8];
        minfreq = 6;
        if istrue(combtrials)
          trialsel = ones(size(predtc));
        else
          trialsel = {valencetc==1, valencetc==2, valencetc==1 & predtc==1, valencetc==1 & predtc==2, valencetc==2 & predtc==1, valencetc==2 & predtc==2};
        end
    end
    
    % prep PAC settings
    ampfreq        = {[minfreq:1:20 24:4:60]  ,70:4:130           ,100:10:200                    };
    phasfreq       = {ampfreq{1}              ,ampfreq{1}         ,ampfreq{1}                    };
    ampfreqtimwin  = {3./ampfreq{1}           ,3./ampfreq{2}      ,ones(size(ampfreq{3})) .* .050};
    phasfreqtimwin = {3./phasfreq{1}          ,3./phasfreq{2}     ,3./phasfreq{3}};
    amptstohfats   = {'no'                    ,'no'               ,'yes'};
    
    % first cut out data
    cfg = [];
    cfg.toilim = toilim;
    datasel = ft_redefinetrial(cfg,data);
    cfg = [];
    cfg.demean  = 'yes';
    cfg.detrend = 'yes';
    cfg.trials  = 'all';
    datasel = ft_preprocessing(cfg,datasel);
    
    % all PACs
    for iset = 1:3
      if ~exist(fn{iset,iperiod},'file')
        
        % get pac
        wplfdata = [];
        for itrlset = 1:numel(trialsel)
          cfg = [];
          cfg.trials          = find(trialsel{itrlset});
          cfg.ampfreq         = ampfreq{iset};
          cfg.ampfreqtimwin   = ampfreqtimwin{iset};
          cfg.phasfreq        = phasfreq{iset};
          cfg.phasfreqtimwin  = phasfreqtimwin{iset};
          cfg.amptstohfats    = amptstohfats{iset};
          cfg.centerampfreq   = centerampfreq;
          cfg.phasfreqnoamp   = 'no';
          cfg.handlespecedges = 'liberal';
          cfg.normalization   = 'withintrials';
          cfg.keeptrials      = 'no';
          wplfdata{itrlset} = rmr_phaseamplitudecoupling(cfg,datasel);
        end
        if numel(wplfdata) == 1
          wplfdata = wplfdata{1};
        end
        
        % save
        save(fn{iset,iperiod},'wplfdata','-v7.3');
        clear wplfdata
      end
      
      % do the two partitions (only if there was no split over trial types
      if istrue(combtrials)
        for ish = 1:2
          currfn = [fn{iset,iperiod}(1:end-4) '_' 'shpart'  num2str(ish) '.mat'];
          if ~exist(currfn,'file')
            
            % get pac
            cfg = [];
            cfg.trials          = ish:2:numel(datasel.trial);
            cfg.ampfreq         = ampfreq{iset};
            cfg.ampfreqtimwin   = ampfreqtimwin{iset};
            cfg.phasfreq        = phasfreq{iset};
            cfg.phasfreqtimwin  = phasfreqtimwin{iset};
            cfg.amptstohfats    = amptstohfats{iset};
            cfg.centerampfreq   = centerampfreq;
            cfg.phasfreqnoamp   = 'no';
            cfg.handlespecedges = 'liberal';
            cfg.normalization   = 'withintrials';
            cfg.keeptrials      = 'no';
            wplfdata = rmr_phaseamplitudecoupling(cfg,datasel);
            
            % save
            save(currfn,'wplfdata','-v7.3');
            clear wplfdata
          end
        end % ish
      end
      
    end % iset
    
  end % iperiod
end














function playground












%%%%%%%%%%%%%%%%%
%%% BETWEEN
% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

% design
centerampfreq = 'yes';
combtrials    = 'no';
for     isubj = [2 3]%1  :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  fn = [];
  fn{1,1} = [info.savepath currsubj '_' 'wplf_4to60'           '_'         'CTI' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials '.mat'];
  fn{2,1} = [info.savepath currsubj '_' 'wplf_4to60to70to130'  '_'         'CTI' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{3,1} = [info.savepath currsubj '_' 'wplf_4to60toHFA'      '_'         'CTI' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{1,2} = [info.savepath currsubj '_' 'wplf_6to60'           '_' 'postface500' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{2,2} = [info.savepath currsubj '_' 'wplf_6to60to70to130'  '_' 'postface500' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{3,2} = [info.savepath currsubj '_' 'wplf_6to60toHFA'      '_' 'postface500' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  
  % for each of the freqs, plot the tfr
  for iperiod = 1:2
    switch iperiod
      case {1} % CTI
        toilim  = [0.2 1.2];
        condname = {'pred','unpred'};
        %condname = {'predunpred',[]};
      case {2} % postface
        toilim  = [1.2 1.8];
        condname = {'fear','neat','fearpred','fearunpred','neutpred','neutunpred'};
        %condname = {'fearneat',[],'fearpredunpred',[],'neutpredunpred',[]};
    end
    
    for iset = [1]
      
      % load that shit up
      load(fn{iset,iperiod})
      for iwplfset = 1:2%:2:numel(wplfdata)
        currwplfdata1 = wplfdata{iwplfset};
        if iset == 3
          currwplfdata2 = wplfdata{iwplfset+1};
          currwplfdata1.wplf = cat(5,currwplfdata1.wplf,currwplfdata2.wplf);
        end
        
        % set some things
        nchan = numel(currwplfdata1.label);
        
        % get layout
        cfg = [];
        cfg.layout = 'ordered';
        tmp = currwplfdata1;
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
        lay.pos    = lay.pos(ind,:);
        lay.width  = lay.width(ind);
        lay.height = lay.height(ind);
        
        
        % plot
        if iperiod == 1
          figname = [currsubj '-' 'CTI'];
        else
          figname = [currsubj '-' 'postface'];
        end
        if numel(currwplfdata1.ampfreq)>2
          figname = [figname '-' 'osc'];
        else
          figname = [figname '-' 'hfa'];
        end
        figname = [figname '-' condname{iwplfset}];
        
        % browse
        cfg = [];
        cfg.layout      = lay;
        cfg.windowname  = figname;
        cfg.plotlabel   = 'yes';
        cfg.plotoutline = 'no';
        rmr_browsewplfdata(cfg,currwplfdata1);
      end
    end % iset
    
  end % iperiod
  
  
end % isubj
%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%
%%% WITHIN
% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

% design
centerampfreq = 'yes';
combtrials    = 'no';
for     isubj = 1  :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  fn = [];
  fn{1,1} = [info.savepath currsubj '_' 'wplf_4to60'           '_'         'CTI' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials '.mat'];
  fn{2,1} = [info.savepath currsubj '_' 'wplf_4to60to70to130'  '_'         'CTI' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{3,1} = [info.savepath currsubj '_' 'wplf_4to60toHFA'      '_'         'CTI' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{1,2} = [info.savepath currsubj '_' 'wplf_6to60'           '_' 'postface500' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{2,2} = [info.savepath currsubj '_' 'wplf_6to60to70to130'  '_' 'postface500' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  fn{3,2} = [info.savepath currsubj '_' 'wplf_6to60toHFA'      '_' 'postface500' '_' 'centerae' centerampfreq '_' 'combtrials' combtrials  '.mat'];
  
  % for each of the freqs, plot the tfr
  for iperiod = 1:2
    switch iperiod
      case {1} % CTI
        toilim  = [0.2 1.2];
        condname = {'pred','unpred'};
      case {2} % postface
        toilim  = [1.2 1.8];
        condname = {'fear','neat','fearpred','fearunpred','neutpred','neutunpred'};
    end
    for iset = [ 3]
      
      % load that shit up
      load(fn{iset,iperiod})
      for iwplfset = 1:numel(wplfdata)
        currwplfdata1 = wplfdata{iwplfset};
        
        % set some things
        nchan = numel(currwplfdata1.label);
        
        % get layout
        cfg = [];
        cfg.layout = 'ordered';
        tmp = currwplfdata1;
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
        lay.pos    = lay.pos(ind,:);
        lay.width  = lay.width(ind);
        lay.height = lay.height(ind);
        
        % flg
        hashfa = size(currwplfdata1.wplf,3)==1;
        
        % plot
        if iperiod == 1
          figname = [currsubj '-' 'CTI'];
        else
          figname = [currsubj '-' 'postface'];
        end
        if numel(currwplfdata1.ampfreq)>2
          figname = [figname '-' 'osc'];
        else
          figname = [figname '-' 'hfa'];
        end
        figname = [figname '-' condname{iwplfset}];
        figure('numbertitle','off','name',figname)
        hold on
        
        % get clim
        clim = [0 0];
        for ichan = 1:nchan
          clim(2) = max(clim(2),max(max(abs(squeeze(currwplfdata1.wplf(ichan,ichan,:,:))))));
        end
        clim = round(clim*100)./100;
        %clim = [0 0.1];
        
        % rearrange wplf
        wplfold   = currwplfdata1.wplf;
        nfreqamp  = size(wplfold,3);
        nfreqphas = size(wplfold,4);
        wplf = zeros(nchan,nfreqamp,nfreqphas);
        for ichan = 1:nchan
          wplf(ichan,:,:) = wplfold(ichan,ichan,:,:);
        end
        
        % fake freqdata
        fakedata = [];
        fakedata.powspctrm = squeeze(abs(wplf));
        if ~hashfa
          fakedata.time   = currwplfdata1.phasfreq;
          fakedata.freq   = currwplfdata1.ampfreq;
          fakedata.dimord = 'chan_freq_time';
        else
          fakedata.freq = currwplfdata1.phasfreq;
          fakedata.dimord = 'chan_freq';
        end
        fakedata.label = currwplfdata1.label;
        
        % plot
        if hashfa
          cfg = [];
          cfg.layout = lay;
          cfg.vlim   = clim;
          cfg.showlabels = 'yes';
          ft_multiplotER(cfg,fakedata);
        else
          cfg = [];
          cfg.layout = lay;
          cfg.zlim   = clim;
          cfg.showlabels = 'yes';
          ft_multiplotTFR(cfg,fakedata);
        end
        
        %       % prep ticks
        %       ampfreq   = wplfdata.ampfreq';
        %       phasfreq  = wplfdata.phasfreq';
        %       nampfreq  = numel(ampfreq);
        %       nphasfreq = numel(phasfreq);
        %       amptickind    = 1:floor(nampfreq/6):nampfreq;
        %       ampticklabel  = num2str(round(ampfreq(amptickind)));
        %       phastickind   = 1:floor(nphasfreq/6):nphasfreq;
        %       phasticklabel = num2str(round(phasfreq(phastickind)));
        %
        %       % draw data
        %       label = wplfdata.label;
        %       nchan = numel(label);
        %       for ichan = 1:nchan
        %         subplot_tight(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan,0.025)
        %         %
        %         C = abs(squeeze(wplfdata.wplf(ichan,ichan,:,:)));
        %         imagesc(C);
        %         axis xy
        %         caxis(clim)
        %         title(wplfdata.label{ichan})
        %         if ichan == nchan
        %           set(gca,'ytick',amptickind,'yticklabel',ampticklabel,'xtick',phastickind,'xticklabel',phasticklabel)
        %           %colorbar('location','eastoutside')
        %         else
        %           set(gca,'ytick',amptickind,'xtick',amptickind,'yticklabel',[],'xticklabel',[])
        %         end
        %       end
        %       subplot_tight(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan+1,0.025)
        %       caxis(clim)
        %       colorbar
        %       axis off
        
      end
    end % iset
    
  end % iperiod
  
  
end % isubj
%%%%%%%%%%%%%%




















