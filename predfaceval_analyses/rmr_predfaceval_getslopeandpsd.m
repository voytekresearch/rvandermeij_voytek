function rmr_predfaceval_getslopeandpsd





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
  fn{1} = [info.savepath currsubj '_' 'slopeandpsd_mtmfft'     '_'         'CTI'   '.mat'];
  fn{2} = [info.savepath currsubj '_' 'slopeandpsd_mtmfft'     '_' 'postface500'   '.mat'];
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
    switch iperiod
      case {1} % CTI
        toilim  = [0.2 1.2];
      case {2} % postface
        toilim  = [1.2 1.8];
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
    
    % get TFR
    cfg = [];
    cfg.channel    = 'all';
    cfg.trials     = 'all';
    cfg.keeptrials = 'yes';
    cfg.pad        = max([ceil(max(cellfun(@numel,datasel.time)) ./ datasel.fsample) 2]);
    cfg.method     = 'mtmfft';
    cfg.output     = 'pow';
    cfg.taper      = 'hanning';
    cfg.foi        = 1:200;
    freqdata = ft_freqanalysis(cfg,datasel);
    
    % save
    save(fn{iperiod},'freqdata','-v7.3');
    clear freqdata
  end % iperiod
  
end













function playground





% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';


for     isubj = 1  :numel(info.subj)
  
  
  % set
  currsubj   = info.subj{isubj};
  fn = [];
  fn{1} = [info.savepath currsubj '_' 'slopeandpsd_mtmfft'     '_'         'CTI'   '.mat'];
  fn{2} = [info.savepath currsubj '_' 'slopeandpsd_mtmfft'     '_' 'postface500'   '.mat'];
  
  % for each of the freqs, plot the tfr
  for iperiod = 1:2
    
    % load that shit up
    load(fn{iperiod})
    % set some things
    nchan  = numel(freqdata.label);
    ntrial = size(freqdata.powspctrm,1);
    
    % get trialdef
    valencetc = freqdata.trialinfo(:,2); % fearful (1) or neutral (2) face trial
    predtc    = freqdata.trialinfo(:,1); % pred (1) or unpred (2) trial
    
    % get layout
    cfg = [];
    cfg.layout = 'ordered';
    tmp = freqdata;
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
    
    % get clim
    ylim = [];
    ylim(1) = min(log10(freqdata.powspctrm(:)));
    ylim(2) = max(log10(freqdata.powspctrm(:)));
    %ylim = ceil(clim);
    %clim = [0 2];
    
%     % plot
%     figure('numbertitle','off','name',figname)
%     for ichan = 1:nchan
%       subplot(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan)
%       currdat = log10(squeeze(freqdata.powspctrm(:,ichan,:)));
%       x = [freqdata.freq freqdata.freq(end:-1:1)];
%       miny = min(currdat,[],1);
%       maxy = max(currdat,[],1);
%       y = [miny maxy(end:-1:1)];
%       patch(log10(x),y,[.7 .7 .7],'edgecolor','none')
%       %plot(freqdata.freq,currdat','color',[.7 .7 .7])
%       hold on
%       plot(log10(freqdata.freq),mean(currdat,1),'cy')
%       set(gca,'ylim',ylim,'xlim',[1 log10(100)],'ytick',[],'xtick',[])
%     end
    
    
    
    % fit between 40 and 60 and plot
    % first fit
    %   chi  = cell(1,2);
    %   offs = cell(1,2);
    %   for ipred = 1:2
    %     currtrialind = find(predtc==ipred);
    %     ntrial = numel(currtrialind);
    %     chi{ipred} = NaN(nchan,ntrial);
    %     offs{ipred} = NaN(nchan,ntrial);
    %     for ichan = 1:nchan
    %       for itrial = 1:ntrial
    %         currdat = squeeze(freqdata.powspctrm(currtrialind(itrial),ichan,find(freqdata.freq==40):find(freqdata.freq==60)));
    %
    %         % get chi
    %         x = log10(40:60);
    %         y = log10(currdat);
    %         out = robustfit(x,y);
    %
    %         % store
    %         chi{ipred}(ichan,itrial)  = out(2);
    %         offs{ipred}(ichan,itrial) = out(1);
    %       end
    %     end
    %   end
    chi  = NaN(nchan,2);
    offs = NaN(nchan,2);
    for ichan = 1:nchan
      pred1 = squeeze(freqdata.powspctrm(predtc==1,ichan,find(freqdata.freq==40):find(freqdata.freq==60)));
      pred2 = squeeze(freqdata.powspctrm(predtc==2,ichan,find(freqdata.freq==40):find(freqdata.freq==60)));
      
      % get chi
      x = log10(40:60);
      y = log10(mean(pred1));
      pred1out = robustfit(x,y);
      y = log10(mean(pred2));
      pred2out = robustfit(x,y);
      
      % store
      chi(ichan,1)  = pred1out(2);
      chi(ichan,2)  = pred2out(2);
      offs(ichan,1) = pred1out(1);
      offs(ichan,2) = pred2out(1);
    end
    
    % plot
    ylim = [floor(min(chi(:))) 0];
    figure('numbertitle','off','name',[figname '-' num2str(ylim(1)) 'to' num2str(ylim(2))])
    for ichan = 1:nchan
      subplot(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan)
      
      % plot
      %     currchimean = [mean(chi{1}(ichan,:)) mean(chi{2}(ichan,:))];
      %     currchistd  = [std(chi{1}(ichan,:)) std(chi{2}(ichan,:))];
      %     barweb(currchimean,currchistd)
      bar([0 chi(ichan,1); 0 chi(ichan,2)]');
      set(gca,'ylim',ylim,'ytick',[],'xtick',[])
      title(freqdata.label{ichan})
    end
    
    
    
  end % iperiod
  
end % isubj
%%%%%%%%%%%%%%












