function rmr_predfaceval_getfourier





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
  fn{1,1} = [info.savepath currsubj '_' 'fourier_4to40hz_3cyc'   '_'              'CTI'   '.mat'];
  fn{1,2} = [info.savepath currsubj '_' 'dataetc_4to40hz_3cyc'   '_'              'CTI'   '.mat'];
  fn{2,1} = [info.savepath currsubj '_' 'fourier_6to40hz_3cyc'   '_'         'postface'   '.mat'];
  fn{2,2} = [info.savepath currsubj '_' 'dataetc_6to40hz_3cyc'   '_'         'postface'   '.mat'];
  if  all(cellfun(@exist,fn(:),repmat({'file'},[numel(fn) 1])))
    continue
  end
  
  % fetch and preprocess data
  cfg = [];
  cfg.demean      = 'yes';
  cfg.detrend     = 'yes';
  cfg.prestim     = 0.3;
  cfg.poststim    = 0.5+0.3;
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
        toilim   = [0 1.2];
        minfreq  = 4;
      case {2} % postface
        toilim   = [1.2 1.9];
        minfreq  = 6;
    end
    
    % first cut out data
    cfg = [];
    cfg.toilim = toilim;
    datasel = ft_redefinetrial(cfg,data);
    cfg = [];
    cfg.demean  = 'yes';
    cfg.detrend = 'yes';
    cfg.trials  = 'all';
    datasel = ft_preprocessing(cfg,datasel);
    
    % all TFRs
    if ~exist(fn{iperiod,1},'file')
      % get fourier
      % get timeoi's for welch tapering
      freqoi    = minfreq:1:40;
      taplength = 3./freqoi;
      % get them
      [dum maxind] = max(cellfun(@numel,datasel.time));
      timeoi = datasel.time{maxind};
      % get fourier
      fourier = roe_fourierwelcheig(datasel,freqoi,timeoi,taplength,'hanning',[]);
      
      % save fourier
      fourfn = fn{iperiod,1};
      save(fourfn,'fourier','-v7.3')
      % save data and stuff
      datasel = rmfield(datasel,'trial');
      datafn = fn{iperiod,2};
      save(datafn,'datasel','freqoi','timeoi')
    end % exist
    
  end % iperiod
end













function playground















% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

for     isubj = 2:5%  :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  
  % check whether already done
  fn = [];
  fn{1,1} = [info.savepath currsubj '_' 'fourier_4to40hz_3cyc'   '_'              'CTI'   '.mat'];
  fn{1,2} = [info.savepath currsubj '_' 'dataetc_4to40hz_3cyc'   '_'              'CTI'   '.mat'];
  fn{2,1} = [info.savepath currsubj '_' 'fourier_6to40hz_3cyc'   '_'         'postface'   '.mat'];
  fn{2,2} = [info.savepath currsubj '_' 'dataetc_6to40hz_3cyc'   '_'         'postface'   '.mat'];
  % skip if needed
  if ~all(cellfun(@exist,fn(:),repmat({'file'},[numel(fn) 1])))
    %continue
  end
  
  % load dataetc to get trialinfo
  load(fn{1,2})
  valencetc = datasel.trialinfo(:,2); % fearful (1) or neutral (2) face trial
  predtc    = datasel.trialinfo(:,1); % pred (1) or unpred (2) trial
  
  % for each of the freqs, plot the tfr
  for iperiod = 1:2
    switch iperiod
      case {1} % CTI
        trialsel = {predtc==1, predtc==2};
        setname  = {'predunpred'};
      case {2} % postface
        trialsel  = {valencetc==1, valencetc==2; valencetc==1 & predtc==1, valencetc==2 & predtc==1; valencetc==1 & predtc==2, valencetc==2 & predtc==2};
        setname  = {'fearneut','fearneutpred','fearneutunpred'};
    end
    
    % load that shit up
    load(fn{iperiod,1})
    load(fn{iperiod,2})
    
    
    % set some things
    nchan = numel(datasel.label);
    
    % do it per set
    for itrialset = 1%:size(trialsel,1)
      ind1 = find(trialsel{itrialset,1});
      ind2 = find(trialsel{itrialset,2});
      
      % get csd
      fourier1 = fourier(:,:,ind1,:);
      fourier1 = roe_fouriernormalize(fourier1,'coh');
      csd1 = zeros(size(fourier1,1),size(fourier1,1),size(fourier1,2));
      for itrial = 1:size(fourier1,3)
        for ifreq = 1:numel(freqoi)
          tmpfour = squeeze(fourier1(:,ifreq,itrial,:));
          tmpfour(:,isnan(tmpfour(1,:))) = [];
          csd1(:,:,ifreq) = csd1(:,:,ifreq) + ((tmpfour * tmpfour')./numel(ind1));
        end
      end
      fourier2 = fourier(:,:,ind2,:);
      fourier2 = roe_fouriernormalize(fourier2,'coh');
      csd2 = zeros(size(fourier2,1),size(fourier2,1),size(fourier2,2));
      for itrial = 1:size(fourier2,3)
        for ifreq = 1:numel(freqoi)
          tmpfour = squeeze(fourier2(:,ifreq,itrial,:));
          tmpfour(:,isnan(tmpfour(1,:))) = [];
          csd2(:,:,ifreq) = csd2(:,:,ifreq) + ((tmpfour * tmpfour')./numel(ind2));
        end
      end
      
      % create fake wplfdata
      wplfdata = [];
      wplfdata.wplf     = cat(5,permute(csd1,[1 2 4 3]),permute(csd2,[1 2 4 3]));
      wplfdata.dimord   = 'ampchan_phaschan_ampfreq_phasfreq';
      wplfdata.ampfreq  = 1;
      wplfdata.phasfreq = freqoi;
      wplfdata.label    = datasel.label;
      
      
      % get layout
      cfg = [];
      cfg.layout = 'ordered';
      tmp = wplfdata;
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
        figname = [currsubj '-' 'CTI' '-'];
      else
        figname = [currsubj '-' 'postface' '-'];
      end
      figname = [figname setname{itrialset}];
      
      % browse
      cfg = [];
      cfg.layout      = lay;
      cfg.windowname  = figname;
      cfg.plotlabel   = 'yes';
      cfg.plotoutline = 'no';
      rmr_browsewplfdata(cfg,wplfdata);
    end
  end % iperiod
  
  
  
end % isubj










































