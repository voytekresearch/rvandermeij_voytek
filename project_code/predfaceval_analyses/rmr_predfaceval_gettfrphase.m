function rmr_predfaceval_gettfrphase





% fetch info
info = rmr_predfaceval_info;

%info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

%
cyclenum = 3;
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
  fn{1} = [info.savepath currsubj '_' 'tfrphase_4to30hz_' num2str(cyclenum) 'cyc'                          '.mat'];
  fn{2} = [info.savepath currsubj '_' 'tfrphase_4to30hz_' num2str(cyclenum) 'cyc'   '_'  'processed'       '.mat'];
  fn{3} = [info.savepath currsubj '_' 'tfrphase_4to30hz_' num2str(cyclenum) 'cyc'   '_'  'processed2'      '.mat'];
  if all(cellfun(@exist,fn(:),repmat({'file'},[numel(fn) 1])))
    continue
  end
  
  if ~exist(fn{1},'file')
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
    
    % get TFR
    cfg = [];
    cfg.channel    = 'all';
    cfg.trials     = 'all';
    cfg.keeptrials = 'yes';
    cfg.pad        = ceil(max(cellfun(@numel,data.time)) ./ data.fsample);
    cfg.method     = 'mtmconvol';
    cfg.output     = 'fourier';
    cfg.taper      = 'hanning';
    cfg.foi        = 4:30;
    cfg.t_ftimwin  = cyclenum./cfg.foi;
    cfg.toi        = 0:(1/250):1.6;
    freqdata = ft_freqanalysis(cfg,data);
    
    % save
    save(fn{1},'freqdata','-v7.3');
    clear data
  else
    load(fn{1})
  end % iset
  
  
  % get trialdef
  valencetc = freqdata.trialinfo(:,2); % fearful (1) or neutral (2) face trial
  predtc    = freqdata.trialinfo(:,1); % pred (1) or unpred (2) trial
  
  % set ns
  nchan = numel(freqdata.label);
  nfreq = numel(freqdata.freq);
  ntime = numel(freqdata.time);
  
  if ~exist(fn{2},'file')
    % compute intertrial phase-locking (coh and plv), and get permutation distr for plv
    % plv
    plvdata = rmfield(freqdata,'fourierspctrm');
    plvdata.dimord = 'chan_freq_time_pred';
    plvdata.plvspctrm = NaN(nchan,nfreq,ntime,2);
    for ipred = 1:2
      trialbool = predtc == ipred;
      for ifreq = 1:nfreq
        currplv = squeeze(freqdata.fourierspctrm(trialbool,:,ifreq,:));
        currplv = abs(squeeze(mean(exp(1i .* angle(currplv)),1)));
        plvdata.plvspctrm(:,ifreq,:,ipred) = currplv;
      end
    end
    
    % coh
    cohdata = rmfield(freqdata,'fourierspctrm');
    cohdata.dimord = 'chan_freq_time_pred';
    cohdata.cohspctrm = NaN(nchan,nfreq,ntime,2);
    for ipred = 1:2
      trialbool = predtc == ipred;
      for ifreq = 1:nfreq
        currcoh = squeeze(freqdata.fourierspctrm(trialbool,:,ifreq,:));
        currcoh = bsxfun(@rdivide, currcoh, sum(abs(currcoh),1));
        currcoh = abs(squeeze(sum(currcoh,1)));
        cohdata.cohspctrm(:,ifreq,:,ipred) = currcoh;
      end
    end
    
    % permute baby!
    coh99spctrm = NaN(nchan,nfreq,ntime,2);
    plv99spctrm = NaN(nchan,nfreq,ntime,2);
    coh95spctrm = NaN(nchan,nfreq,ntime,2);
    plv95spctrm = NaN(nchan,nfreq,ntime,2);
    nperm   = 1000;
    permlim = round(1 .* 250); % 'sampling' is now 250Hz
    for ipred = 1:2
      trialbool = predtc == ipred;
      currfour  = squeeze(freqdata.fourierspctrm(trialbool,:,:,:));
      ntrial    = sum(trialbool);
      for ichan = 1:nchan
        disp(['obtaining permutation for channel ' num2str(ichan) '/' num2str(nchan) ' pred' num2str(ipred)])
        tmpcoh99spctrm = NaN(nfreq,ntime);
        tmpplv99spctrm = NaN(nfreq,ntime);
        tmpcoh95spctrm = NaN(nfreq,ntime);
        tmpplv95spctrm = NaN(nfreq,ntime);
        
        for ifreq = 1:nfreq
          currdat     = permute(currfour(:,ichan,ifreq,:),[1 4 2 3]);
          nonnanbool   = ~sum(isnan(currdat),1);
          currdat     = currdat(:,nonnanbool);
          nnnantime   = sum(nonnanbool);
          permdistplv = zeros(nperm,nnnantime);
          permdistcoh = zeros(nperm,nnnantime);
          for iperm = 1:nperm
            permdat = currdat;
            for itrial = 1:ntrial
              currshift = round(rand(1) .* permlim);
              permdat(itrial,:) = circshift(permdat(itrial,:),currshift,2);
            end
            % calc coh and plv
            permdistcoh(iperm,:) = abs(sum(bsxfun(@rdivide,permdat,sum(abs(permdat),1)),1));
            %permdistplv(iperm,:) = abs(mean(exp(1i.*angle(permdat))));
            permdistplv(iperm,:) = abs(mean(permdat ./ abs(permdat),1));
          end
          % get 99th percentile
          permdistcoh = sort(permdistcoh,1);
          permdistplv = sort(permdistplv,1);
          % insert
          tmpcoh99spctrm(ifreq,nonnanbool) = permdistcoh(min([round(nperm .* .99)+1 nperm]),:);
          tmpcoh95spctrm(ifreq,nonnanbool) = permdistcoh(min([round(nperm .* .99)+1 nperm]),:);
          tmpplv99spctrm(ifreq,nonnanbool) = permdistplv(min([round(nperm .* .99)+1 nperm]),:);
          tmpplv95spctrm(ifreq,nonnanbool) = permdistplv(min([round(nperm .* .99)+1 nperm]),:);
        end
        coh99spctrm(ichan,:,:,ipred) = tmpcoh99spctrm;
        plv99spctrm(ichan,:,:,ipred) = tmpcoh95spctrm;
        coh95spctrm(ichan,:,:,ipred) = tmpplv99spctrm;
        plv95spctrm(ichan,:,:,ipred) = tmpplv95spctrm;
      end
    end
    
    % insert
    cohdata.coh99spctrm = coh99spctrm;
    cohdata.coh95spctrm = coh95spctrm;
    plvdata.plv99spctrm = plv99spctrm;
    plvdata.plv95spctrm = plv95spctrm;
    
    % get wwspctrm's, put in wwdata
    wwdata = rmfield(freqdata,'fourierspctrm');
    wwdata.dimord        = 'chan_freq_time';
    wwdata.wwspctrm      = NaN(nchan,nfreq,ntime);
    wwdata.wwnpspctrm    = NaN(nchan,nfreq,ntime);
    wwdata.assviolspctrm = false(nchan,nfreq,ntime);
    for ifreq = 1:nfreq
      disp(['obtaining ww and wwnp freq ' num2str(ifreq) '/' num2str(nfreq)])
      phasdat = angle(squeeze(freqdata.fourierspctrm(:,:,ifreq,:)));
      %ampdat  = abs(squeeze(freqdata.fourierspctrm(:,:,ifreq,:)));
      for ichan = 1:nchan
        for itime = 1:ntime
          currphas = phasdat(:,ichan,itime);
          if ~any(isnan(currphas))
            [dum,table,assviol] = roe_circ_wwtest(currphas,valencetc);
            wwdata.wwspctrm(ichan,ifreq,itime)      = sqrt(table{2,5});
            wwdata.assviolspctrm(ichan,ifreq,itime) = assviol;
            [dum,dum,stat] = circ_cmtest(currphas,valencetc);
            wwdata.wwnpspctrm(ichan,ifreq,itime)    = stat;
          end
        end
      end
    end
    
    % compute rtestdata
    rtestdata = rmfield(freqdata,'fourierspctrm');
    rtestdata.dimord = 'chan_freq_time_pred';
    rtestdata.rtestzspctrm = NaN(nchan,nfreq,ntime,2);
    rtestdata.rtestpspctrm = NaN(nchan,nfreq,ntime,2);
    for ipred = 1:2
      trialbool = predtc == ipred;
      for ichan = 1:nchan
        for ifreq = 1:nfreq
          for itime = 1:ntime
            currangle = angle(squeeze(freqdata.fourierspctrm(trialbool,ichan,ifreq,itime)));
            if any(isnan(currangle))
              continue
            end
            [pval zval] = circ_rtest(currangle);
            rtestdata.rtestzspctrm(ichan,ifreq,itime,ipred) = zval;
            rtestdata.rtestpspctrm(ichan,ifreq,itime,ipred) = pval;
          end
        end
      end
    end
    
    % save that shit
    save(fn{2},'cohdata','plvdata','wwdata','rtestdata','-v7.3');
    clear freqdata cohdata plvdata wwdata rtestdata
  else
    load(fn{2})
  end % iset
  
  
  if ~exist(fn{3},'file')
    % compute plv and cohdata condition difference using permutation testing based on the Z statistic (Brillinger 1981 and Amjad et al 1997)
    ntrialpred   = sum(predtc==1);
    ntrialunpred = sum(predtc==2);
    ntrial       = ntrialpred+ntrialunpred;
    % prep coh
    cohdiffdata = rmfield(cohdata,{'coh99spctrm','coh95spctrm','cohspctrm'});
    cohdiffdata.dimord         = 'chan_freq_time';
    cohdiffdata.cohdiffspctrm  = cohdata.cohspctrm(:,:,:,1)-cohdata.cohspctrm(:,:,:,2);
    cohdiffdata.cohdiffZspctrm = ((atanh(cohdata.cohspctrm(:,:,:,1))-(1/(ntrialpred-1))) -  (atanh(cohdata.cohspctrm(:,:,:,2))-(1/(ntrialunpred-1)))) ./ sqrt((1/(ntrialpred-1))+(1/(ntrialunpred-1)));
    % prep plv
    plvdiffdata = rmfield(plvdata,{'plv99spctrm','plv95spctrm','plvspctrm'});
    plvdiffdata.dimord         = 'chan_freq_time';
    plvdiffdata.plvdiffspctrm  = plvdata.plvspctrm(:,:,:,1)-plvdata.plvspctrm(:,:,:,2);
    plvdiffdata.plvdiffZspctrm = ((atanh(plvdata.plvspctrm(:,:,:,1))-(1/(ntrialpred-1))) -  (atanh(plvdata.plvspctrm(:,:,:,2))-(1/(ntrialunpred-1)))) ./ sqrt((1/(ntrialpred-1))+(1/(ntrialunpred-1)));
    
    % get permutation critvals
    cohZpercspctrm = NaN(nchan,nfreq,ntime);
    plvZpercspctrm = NaN(nchan,nfreq,ntime);
    nperm   = 1000;
    for ichan = 1:nchan
      disp(['obtaining permutation for channel ' num2str(ichan) '/' num2str(nchan)])
      tmpcohZpercspctrm = NaN(nfreq,ntime);
      tmpplvZpercspctrm = NaN(nfreq,ntime);
      
      for ifreq = 1:nfreq
        currdat     = permute(freqdata.fourierspctrm(:,ichan,ifreq,:),[1 4 2 3]);
        nonnanbool   = ~sum(isnan(currdat),1);
        currdat     = currdat(:,nonnanbool);
        nnnantime   = sum(nonnanbool);
        permdistZplv = zeros(nperm,nnnantime);
        permdistZcoh = zeros(nperm,nnnantime);
        for iperm = 1:nperm
          permpredind = false(1,ntrial);
          permpredind(randperm(ntrial,ntrialpred)) = true;
          currpreddat   = currdat(permpredind,:);
          currunpreddat = currdat(~permpredind,:);
          % calc coh and plv
          currcohpred   = abs(sum(bsxfun(@rdivide,currpreddat,sum(abs(currpreddat),1)),1));
          currcohunpred = abs(sum(bsxfun(@rdivide,currunpreddat,sum(abs(currunpreddat),1)),1));
          currplvpred   = abs(mean(currpreddat ./ abs(currpreddat),1));
          currplvunpred = abs(mean(currunpreddat ./ abs(currunpreddat),1));
          % get Zs
          permdistZcoh(iperm,:) = ((atanh(currcohpred)-(1/(ntrialpred-1))) -  (atanh(currcohunpred)-(1/(ntrialunpred-1)))) ./ sqrt((1/(ntrialpred-1))+(1/(ntrialunpred-1)));
          permdistZplv(iperm,:) = ((atanh(currplvpred)-(1/(ntrialpred-1))) -  (atanh(currplvunpred)-(1/(ntrialunpred-1)))) ./ sqrt((1/(ntrialpred-1))+(1/(ntrialunpred-1)));
        end
        % compute percentile the data is in
        % coh
        permdistZcoh = sort(permdistZcoh,1);
        permindZcoh = bsxfun(@minus,permdistZcoh,permute(cohdiffdata.cohdiffZspctrm(ichan,ifreq,nonnanbool),[1 3 2]));
        permindZcoh(permindZcoh<0) = inf;
        [dum, permindZcoh] = min(permindZcoh,[],1);
        % plv
        permdistZplv = sort(permdistZplv,1);
        permindZplv = bsxfun(@minus,permdistZplv,permute(plvdiffdata.plvdiffZspctrm(ichan,ifreq,nonnanbool),[1 3 2]));
        permindZplv(permindZplv<0) = inf;
        [dum, permindZplv] = min(permindZplv,[],1);
        
        % insert
        tmpcohZpercspctrm(ifreq,nonnanbool) = permindZcoh ./ nperm;
        tmpplvZpercspctrm(ifreq,nonnanbool) = permindZplv ./ nperm;
      end
      cohZpercspctrm(ichan,:,:) = tmpcohZpercspctrm;
      plvZpercspctrm(ichan,:,:) = tmpplvZpercspctrm;
    end
    
    % insert
    cohdiffdata.cohdiffZpercspctrm = cohZpercspctrm;
    plvdiffdata.plvdiffZpercspctrm = plvZpercspctrm;
    
    % save that shit
    save(fn{3},'cohdiffdata','plvdiffdata','-v7.3');
    clear freqdata cohdata plvdata wwdata rtestdata cohdiffdata plvdiffdata
  end
  
  
  
end















function playground





% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

cyclenum = 3;
for     isubj = 1  :numel(info.subj)
  
  % set
  currsubj   = info.subj{isubj};
  fn = [];
  fn{1} = [info.savepath currsubj '_' 'tfrphase_4to30hz_' num2str(cyclenum) 'cyc'                          '.mat'];
  fn{2} = [info.savepath currsubj '_' 'tfrphase_4to30hz_' num2str(cyclenum) 'cyc'   '_'  'processed'       '.mat'];
  fn{3} = [info.savepath currsubj '_' 'tfrphase_4to30hz_' num2str(cyclenum) 'cyc'   '_'  'processed2'      '.mat'];
  
  
  % load that shit up
  %load(fn{1})
  load(fn{2})
  load(fn{3})
  
  
  
  % get trialdef
  valencetc = plvdata.trialinfo(:,2); % fearful (1) or neutral (2) face trial
  predtc    = plvdata.trialinfo(:,1); % pred (1) or unpred (2) trial
  predlabel = {'pred','unpred'};
  
  % set some things
  nchan  = numel(plvdata.label);
  nfreq  = size(plvdata.plvspctrm,2);
  ntime  = size(plvdata.plvspctrm,3);
  
  % get layout
  cfg = [];
  cfg.layout = 'ordered';
  tmp = wwdata;
  tmp.label = tmp.label;
  lay = ft_prepare_layout(cfg,tmp);
  
  % reorder channels in lay
  ind = [];
  type = {'RAM','LAM','AM','ROF','LOF','OF','RIN','LIN','IN'};
  for itype = 1:numel(type)
    ind = [ind; find(strncmp(lay.label,type{itype},numel(type{itype})))];
  end
  chanind = [ind; setdiff(1:nchan,ind)'];
  lay.label  = lay.label(chanind);
  %lay.pos    = lay.pos(ind,:);
  lay.width  = lay.width(chanind);
  lay.height = lay.height(chanind);
  
  
  
  
  %
  %   % plot
  %   for ipred = 1:2
  %     tmpdata =                        cohdata;
  %     tmpdata.interestfield1 = tmpdata.cohspctrm;
  %     tmpdata.interestfield2 = tmpdata.coh95spctrm;
  %     %
  %     tmpdata.interestfield1 = tmpdata.interestfield1(:,:,:,ipred);
  %     tmpdata.interestfield1(tmpdata.interestfield1<tmpdata.interestfield2(:,:,:,ipred)) = NaN;
  %     tmpdata.interestfield1(isnan(tmpdata.interestfield2(:,:,:,ipred))) = NaN;
  %     tmpdata.dimord    = 'chan_freq_time';
  %     figure('numbertitle','off','name',[currsubj '-' predlabel{ipred}])
  %     for ichan = 1:nchan
  %       subplot(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan)
  %       cfg = [];
  %       cfg.channel     = chanind(ichan);
  %       cfg.parameter   = 'interestfield1';
  %       cfg.colorbar    = 'no';
  %       cfg.ylim        = [4 15];
  %       cfg.zlim        = [0 1];
  %       cfg.interactive = 'no';
  %       ft_singleplotTFR(cfg,tmpdata);
  %       axis off
  %       % lines
  %       ylim = get(gca,'ylim');
  %       xlim = get(gca,'xlim');
  %       line([0 0],ylim,'color','w','linewidth',2)
  %       line([.2 .2],ylim,'color','w','linewidth',2)
  %       line([1.2 1.2],ylim,'color','w','linewidth',2)
  %       % axes
  %       line(xlim,[ylim(1) ylim(1)],'color','k')
  %       line(xlim,[ylim(2) ylim(2)],'color','k')
  %       line([xlim(1) xlim(1)],ylim,'color','k')
  %       line([xlim(2) xlim(2)],ylim,'color','k')
  %     end
  %   end
  %
  %
  %
  
  %   % plot
  %   for ipred = 1:2
  %     tmpdata =                        rtestdata;
  %     tmpdata.interestfield1 = tmpdata.rtestzspctrm(:,:,:,ipred);
  %     tmpdata.dimord    = 'chan_freq_time';
  %     figure('numbertitle','off','name',[currsubj '-' 'R' predlabel{ipred}])
  %     for ichan = 1:nchan
  %       subplot(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan)
  %       cfg = [];
  %       cfg.channel     = chanind(ichan);
  %       cfg.parameter   = 'interestfield1';
  %       cfg.colorbar    = 'no';
  %       cfg.ylim        = [4 15];
  %       cfg.zlim        = [0 5];
  %       cfg.interactive = 'no';
  %       ft_singleplotTFR(cfg,tmpdata);
  %       axis off
  %       % lines
  %       ylim = get(gca,'ylim');
  %       xlim = get(gca,'xlim');
  %       line([0 0],ylim,'color','w','linewidth',2)
  %       line([.2 .2],ylim,'color','w','linewidth',2)
  %       line([1.2 1.2],ylim,'color','w','linewidth',2)
  %       % axes
  %       line(xlim,[ylim(1) ylim(1)],'color','k')
  %       line(xlim,[ylim(2) ylim(2)],'color','k')
  %       line([xlim(1) xlim(1)],ylim,'color','k')
  %       line([xlim(2) xlim(2)],ylim,'color','k')
  %     end
  %   end
  
  
  % plot ZDATA
  tmpdata = cohdiffdata;
  tmpdata.interestfield1 = tmpdata.cohdiffZspctrm;
  tmpdata.interestfield2 = tmpdata.cohdiffZpercspctrm;
  %
  tmpdata.interestfield1 = tmpdata.interestfield1(:,:,:);
  tmpdata.interestfield1(tmpdata.interestfield2<.95) = NaN;
  tmpdata.interestfield1(isnan(tmpdata.interestfield2(:,:,:))) = NaN;
  tmpdata.dimord    = 'chan_freq_time';
  figure('numbertitle','off','name',[currsubj '-z'])
  for ichan = 1:nchan
    subplot(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan)
    cfg = [];
    cfg.channel     = chanind(ichan);
    cfg.parameter   = 'interestfield1';
    cfg.colorbar    = 'no';
    cfg.ylim        = [4 15];
    cfg.zlim        = [-2 2];
    cfg.interactive = 'no';
    ft_singleplotTFR(cfg,tmpdata);
    axis off
    % lines
    ylim = get(gca,'ylim');
    xlim = get(gca,'xlim');
    %line([0 0],ylim,'color','w','linewidth',2)
    line([.2 .2],ylim,'color','w','linewidth',2)
    line([1.2 1.2],ylim,'color','w','linewidth',2)
    % axes
    line(xlim,[ylim(1) ylim(1)],'color','k')
    line(xlim,[ylim(2) ylim(2)],'color','k')
    line([xlim(1) xlim(1)],ylim,'color','k')
    line([xlim(2) xlim(2)],ylim,'color','k')
  end
  %   figure
  %   cfg = [];
  %   cfg.channel     = chanind(ichan);
  %   cfg.parameter   = 'interestfield1';
  %   cfg.colorbar    = 'no';
  %   cfg.ylim        = [4 15];
  %   cfg.zlim        = [-2 2];
  %   cfg.interactive = 'no';
  %   ft_singleplotTFR(cfg,tmpdata);
  %   % lines
  %   ylim = get(gca,'ylim');
  %   xlim = get(gca,'xlim');
  %   %line([0 0],ylim,'color','w','linewidth',2)
  %   line([.2 .2],ylim,'color','k','linewidth',2)
  %   line([1.2 1.2],ylim,'color','k','linewidth',2)
  %   % axes
  %   line(xlim,[ylim(1) ylim(1)],'color','k')
  %   line(xlim,[ylim(2) ylim(2)],'color','k')
  %   line([xlim(1) xlim(1)],ylim,'color','k')
  %   line([xlim(2) xlim(2)],ylim,'color','k')
  %   ylabel('frequency (Hz)')
  %   xlabel('time (s)')
  %   set(gca,'xlim',xlim,'xtick',[0 .2 1.2])
  %   set(gca,'ylim',ylim,'ytick',4:2:15)
  %   % axes
  %   colorbar
  %   caxis([-2 2])
  
  
  %   % plot
  %   tmpdata = wwdata;
  %   tmpdata.interestfield1 = tmpdata.wwnpspctrm;
  %   tmpdata.dimord    = 'chan_freq_time';
  %   figure('numbertitle','off','name',[currsubj '-ww'])
  %   for ichan = 1:nchan
  %     subplot(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan)
  %     cfg = [];
  %     cfg.channel     = chanind(ichan);
  %     cfg.parameter   = 'interestfield1';
  %     cfg.colorbar    = 'no';
  %     cfg.ylim        = [4 15];
  %     cfg.zlim        = [0 10];
  %     cfg.interactive = 'no';
  %     ft_singleplotTFR(cfg,tmpdata);
  %     axis off
  %     % lines
  %     ylim = get(gca,'ylim');
  %     xlim = get(gca,'xlim');
  %     line([0 0],ylim,'color','w','linewidth',2)
  %     line([.2 .2],ylim,'color','w','linewidth',2)
  %     line([1.2 1.2],ylim,'color','w','linewidth',2)
  %     % axes
  %     line(xlim,[ylim(1) ylim(1)],'color','k')
  %     line(xlim,[ylim(2) ylim(2)],'color','k')
  %     line([xlim(1) xlim(1)],ylim,'color','k')
  %     line([xlim(2) xlim(2)],ylim,'color','k')
  %   end
  %
  
  
  
  
  
  %   % set lims
  %   freqind = 5; % 8Hz
  %   timelim = [0.8 1.1];
  %   timeind = nearest(freqdata.time,timelim(1)):nearest(freqdata.time,timelim(2));
  %   %timeind = nearest(freqdata.time,0.9);
  %   % plot
  %   for ipred = 1:2
  %     figure('numbertitle','off','name',[currsubj predlabel{ipred}])
  %     for ichan = 1:nchan
  %       subplot(ceil(sqrt(nchan)),ceil(sqrt(nchan)),ichan)
  %       % get
  %       %plotdat = angle(mean(squeeze(freqdata.fourierspctrm(predtc==ipred,chanind(ichan),freqind,timeind)),2));
  %       plotdat = angle(mean(exp(1i.*angle(squeeze(freqdata.fourierspctrm(predtc==ipred,chanind(ichan),freqind,timeind)))),2));
  %       roe_rose(plotdat);
  %       hold on
  %       %     h = roe_compass(mean(plotdat),'r');
  %       %     set(h,'linewidth',2)
  %       title(wwdata.label{chanind(ichan)})
  %     end
  %   end
  
  
  
  
end % isubj
%%%%%%%%%%%%%%












































%
%
% % fetch info
% info = rmr_predfaceval_info;
%
% % info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
% % info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';
%
% for     isubj = 1  :numel(info.subj)
%   % set
%   currsubj   = info.subj{isubj};
%
%   % for artifacts been detected
%   if ~info.(currsubj).trlartfctflg
%     disp(['artifacts not yet detected for ' currsubj])
%     continue
%   end
%
%
%
%   % check whether already done
%   fn = [];
%   fn{1} = [info.savepath currsubj '_' 'tfrphase_4to30hz_3cyc'                          '.mat'];
%   fn{2} = [info.savepath currsubj '_' 'tfrphase_4to30hz_3cyc'   '_'  'processed'       '.mat'];
%
%
%   load(fn{1})
%   %load(fn{2})
%
%   if exist(fn{2},'file')
%     % get trialdef
%     valencetc = freqdata.trialinfo(:,2); % fearful (1) or neutral (2) face trial
%     predtc    = freqdata.trialinfo(:,1); % pred (1) or unpred (2) trial
%
%     % set ns
%     nchan = numel(freqdata.label);
%     nfreq = numel(freqdata.freq);
%     ntime = numel(freqdata.time);
%
%     % compute rtestdata
%     rtestdata = rmfield(freqdata,'fourierspctrm');
%     rtestdata.dimord = 'chan_freq_time_pred';
%     rtestdata.rtestzspctrm = NaN(nchan,nfreq,ntime,2);
%     rtestdata.rtestpspctrm = NaN(nchan,nfreq,ntime,2);
%     for ipred = 1:2
%       trialbool = predtc == ipred;
%       for ichan = 1:nchan
%         for ifreq = 1:nfreq
%           for itime = 1:ntime
%             currangle = angle(squeeze(freqdata.fourierspctrm(trialbool,ichan,ifreq,itime)));
%             if any(isnan(currangle))
%               continue
%             end
%             [pval zval] = circ_rtest(currangle);
%             rtestdata.rtestzspctrm(ichan,ifreq,itime,ipred) = zval;
%             rtestdata.rtestpspctrm(ichan,ifreq,itime,ipred) = pval;
%           end
%         end
%       end
%     end
%
%     % save that shit
%     load(fn{2});
%     save(fn{2},'cohdata','plvdata','wwdata','rtestdata','-v7.3');
%     clear freqdata cohdata plvdata wwdata rtestdata
%   end
% end
%






