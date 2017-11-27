function rmr_abcaudecog_getpacdecomp



% get details
info = rmr_abcaudecog_info;
%info.datapath = '/Users/roemer/Work/Data/tmpdata/BidetCaulet_ECoG/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/abcaudecog/';


% sub select for amua differences between conditions of some sort
%info.subj = {'GP15','GP22','GP28','GP35'};

% nway
nwayalg     = 'parafac';
nwaynmethod = 'degeneracy'; % 'corcondiag'  'splithalf'   'degeneracy'
nwayadd     = 'rnd100_conv1e-6'; %
%
%comptrlload = 0;

% set
centerampfreq = 'yes';

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp('*******')
  disp('*******')
  disp('*******')
  disp('*******')
  disp('*******')
  disp('*******')
  disp(upper(['*******   working on ' currsubj ' - ' nwayalg '_' nwaynmethod '_' nwayadd]))
  disp('*******')
  disp('*******')
  disp('*******')
  disp('*******')
  disp('*******')
  disp('*******')
  
  %
  datfn = [];
  datfn{1} = [info.savepath currsubj '_' 'wplf_9to60'           '_' 'centerae' centerampfreq  '.mat'];
  datfn{2} = [info.savepath currsubj '_' 'wplf_9to60to70to130'  '_' 'centerae' centerampfreq  '.mat'];
  datfn{3} = [info.savepath currsubj '_' 'wplf_9to60toHFA'      '_' 'centerae' centerampfreq  '.mat'];
  %
  nwayfn = cell(numel(datfn),1);
  for inway = 1:numel(datfn)
    nwayfn{inway} = [datfn{inway}(1:end-4) '_' nwayalg '_' nwaynmethod '_' nwayadd '.mat'];
  end
  %
  for idataset = 1:numel(datfn)
    if ~exist(datfn{idataset},'file')
      continue
    end
    if ~exist(nwayfn{idataset},'file')
      
      % deal with saving the data in a file if necessary
      if any(strcmp(nwaynmethod,{'corcondiag','degeneracy'}))
        load(datfn{idataset})
        wplf = squeeze(wplfdata.wplf);
        if false;%ndims(wplf) == 4
          datonlyfn = [datfn{idataset}(1:end-4) '_' 'wplfsonly.mat'];
          save(datonlyfn,'wplf','-v7.3');
          wplfdata.wplf = datonlyfn;
        else
          wplfdata.wplf = wplf;
        end
      elseif strcmp(nwaynmethod,'splithalf')
        sh1fn = [datfn{idataset}(1:end-4) '_shpart1' '.mat'];
        sh2fn = [datfn{idataset}(1:end-4) '_shpart2' '.mat'];
        sh1   = load(sh1fn);
        sh2   = load(sh2fn);
        load(datfn{idataset})
        wplf    = squeeze(wplfdata.wplf);
        wplfsh1 = squeeze(sh1.wplfdata.wplf);
        wplfsh2 = squeeze(sh2.wplfdata.wplf);
        if false;%ndims(wplf) == 4
          datonlyfn    = [datfn{idataset}(1:end-4) '_' 'wplfsonly.mat'];
          sh1datonlyfn = [sh1fn(1:end-4) '_' 'wplfsonly.mat'];
          sh2datonlyfn = [sh2fn(1:end-4) '_' 'wplfsonly.mat'];
          save(datonlyfn,'wplf','-v7.3');
          save(sh1datonlyfn,'wplfsh1','-v7.3');
          save(sh2datonlyfn,'wplfsh2','-v7.3');
          wplfdata.wplf = datonlyfn;
          wplfdata.wplfpart{1} = sh1datonlyfn;
          wplfdata.wplfpart{2} = sh2datonlyfn;
        else
          wplfdata.wplf = wplf;
          wplfdata.wplfpart{1} = wplfsh1;
          wplfdata.wplfpart{2} = wplfsh2;
        end
      else
        error('woops')
      end
      
      
      % nway settings
      cfg = [];
      cfg.model              = nwayalg;         % the model for extracting rhythmic components
      cfg.datparam           = 'wplf';      % the field containing our Fourier coefficients
      cfg.ncompestshdatparam = 'wplfpart';
      cfg.ncompest           = nwaynmethod;        % estimate number of components
      cfg.ncompeststart      = 5;                  % start from 5
      cfg.ncompeststep       = 2;                  % increase number in steps of 2
      cfg.ncompestend        = 25;                 % extract no more than 50
      cfg.ncompestcorconval  = 0.7;
      cfg.numiter            = 1000;               % max number of iteration
      cfg.convcrit           = 1e-6;               % stop criterion of algorithm: minimum relative difference in fit between iterations
      cfg.randstart          = 100;                 % number of random starts
      cfg.ncompestrandstart  = 100;                 % number of random starts for the split-half procedure
      cfg.outputfile         = nwayfn{idataset};
      %
      if ndims(wplf) == 4
        cfg.complexdims        = [1 1 0 0];
        cfg.ncompestshcritval  = [.5 .5 .5 .5];     % splithalf criterion for each of the parameters
      else
        cfg.complexdims        = [1 1 0];
        cfg.ncompestshcritval  = [.5 .5 .5];     % splithalf criterion for each of the parameters
      end
      if false;%ndims(wplf) == 4
        % qsub options
        cfg.distcomp.system          = 'torque';
        cfg.distcomp.timreq          = 60*60*24*.5; % 0.5 days
        %cfg.distcomp.inputsaveprefix = '/home/roevdmei/qsubtemp/nwaytemp_';
        cfg.distcomp.matlabcmd       = '/opt/matlab/2014b/bin/matlab';
        cfg.distcomp.torquequeue     = 'hotel';
      end
      %cfg.distcomp.system          = 'matlabpct';
      %
      nwaycomp = nd_nwaydecomposition(cfg,wplfdata);
    end
    
    
    
    
    
    
    %%% OBTAIN TRIAL LOADINGS
    for comptrlload = 1
      if ~exist([nwayfn{idataset}(1:end-4) '_' 'trialloadings' '_' 'comp' num2str(comptrlload) '.mat'],'file')
        load(datfn{idataset})
        load(nwayfn{idataset})
        
        % fetch and preprocess data
        cfg = [];
        cfg.demean      = 'yes';
        cfg.prestim     = 0;
        cfg.poststim    = 0;
        cfg.reref       = 'yes';
        cfg.refchannel  = 'all';
        data = rmr_abcaudecog_fetchpreprocessdata(cfg,currsubj,info);
        allstdcr   = data.trialinfo(:,2) == 1 & data.trialinfo(:,4) == 4; % standard,  correctly rejected
        
        % cut out data
        cfg = [];
        cfg.toilim = [0 0.4];
        datasel = ft_redefinetrial(cfg,data);
        cfg = [];
        cfg.demean  = 'yes';
        cfg.detrend = 'yes';
        cfg.trials  = find(allstdcr & ((cellfun(@numel,data.time).'./datasel.fsample)>0.380)); % lame
        datasel = ft_preprocessing(cfg,datasel);
        
        % fetch pac cfg
        wplfcfg = wplfdata.cfg;
        
        % get trial loadings
        nwaycomp = rmr_trialloadingpacparafac(nwaycomp,datasel,wplfcfg,comptrlload);
        
        % save
        save([nwayfn{idataset}(1:end-4) '_' 'trialloadings' '_' 'comp' num2str(comptrlload) '.mat'],'nwaycomp','-v7.3');
        clear data datasel nwaycomp wplfdata
      end
    end
    %%%%%
    
    
    
    
    
  end % idataset
end















































function playground





% get details
info = rmr_abcaudecog_info;
%info.datapath = '/Users/roemer/Work/Data/tmpdata/BidetCaulet_ECoG/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/abcaudecog/';


% sub select for amua differences between conditions of some sort
%info.subj = {'GP15','GP22','GP28','GP35'};

% nway
nwayalg     = 'parafac';
nwaynmethod = 'corcondiag'; % 'corcondiag'  'splithalf'
nwayadd     = 'rnd100_conv1e-6'; %

% set
centerampfreq = 'yes';

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1  :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  
  %
  datfn = [];
  datfn{1} = [info.savepath currsubj '_' 'wplf_9to60'           '_' 'centerae' centerampfreq  '.mat'];
  datfn{2} = [info.savepath currsubj '_' 'wplf_9to60to70to130'  '_' 'centerae' centerampfreq  '.mat'];
  datfn{3} = [info.savepath currsubj '_' 'wplf_9to60toHFA'      '_' 'centerae' centerampfreq  '.mat'];
  %
  for idataset = 3%:numel(datfn)
    for icomptrlload = 1
      
      nwayfn = [datfn{idataset}(1:end-4) '_' nwayalg '_' nwaynmethod '_' nwayadd '_' 'trialloadings'  '_' 'comp' num2str(icomptrlload) '.mat'];
      if ~exist(nwayfn,'file')
        continue
      end
      
      load(nwayfn)
      switch idataset
        case {1,2}
          namesuff = 'osc';
        case {3}
          namesuff = 'hfa';
      end
      
      
      
      % plot decomp stats
      if false
        cfg = [];
        cfg.savefig    = 'no';
        cfg.figvisible = 'on';
        roe_nwaystatplot(cfg,nwaycomp)
      end
      
      % set flags
      hashfa = numel(nwaycomp.ampfreq)==1;
      if (hashfa && numel(nwaycomp.comp{1})==4) || (~hashfa && numel(nwaycomp.comp{1})==5)
        hasrpt = true;
      else
        hasrpt = false;
      end
      
      
      % set basics
      comp   = nwaycomp.comp;
      ncomp  = numel(comp);
      nchan  = size(comp{1}{1},1);
      nfreq  = size(comp{1}{3},1);
      if hasrpt
        if hashfa
          ntrial = size(comp{1}{4},1);
        else
          ntrial = size(comp{1}{5},1);
        end
      end
      
      % (all are already correctly rejected, monaural and standard)
      % get attstd/ignstd (monoaural att only)
      % get binaural/monaural attention (all monaural stimuli)
      % get ipsi/contra   (stim, monaural)
      % get lo/hi pitch
      att = (nwaycomp.trialinfo(:,1) == 1 & nwaycomp.trialinfo(:,6) == 1) | (nwaycomp.trialinfo(:,1) == 2 & nwaycomp.trialinfo(:,6) == 2);
      ign = (nwaycomp.trialinfo(:,1) == 1 & nwaycomp.trialinfo(:,6) == 2) | (nwaycomp.trialinfo(:,1) == 2 & nwaycomp.trialinfo(:,6) == 1);
      mon = nwaycomp.trialinfo(:,1) ~= 3;
      bin = nwaycomp.trialinfo(:,1) == 3;
      lo  = nwaycomp.trialinfo(:,7) == 1;
      hi  = nwaycomp.trialinfo(:,7) == 2;
      switch info.(currsubj).leftrightelec
        case 'right'
          ipsi   = nwaycomp.trialinfo(:,6) == 2;
          contra = nwaycomp.trialinfo(:,6) == 1;
        case 'left'
          ipsi   = nwaycomp.trialinfo(:,6) == 1;
          contra = nwaycomp.trialinfo(:,6) == 2;
        otherwise
          error('woopsie')
      end
      
      
      %   1)  4th column of trl: condition (attL = 1, attR = 2, attB = 3)
      %   2)  5th column of trl: standard/deviant (std = 1, dev = 2, 3rd cond = NaN)
      %   3)  6th column of trl: target/no-target (1/0)
      %   4)  7th column of trl: hit/miss/FA/CR (1/2/3/4)
      %   5)  8th column of trl: RT
      %   6)  9th column of trl: left/right/bin stim (left = 1, right = 2, binaural = 3)
      %   7) 10th column of trl: low/high pitch (1/2)
      %   8) 11th column of trl: session/block number
      %   9) 12th column of trl: trial number in order of appearance
      %  10) 13th column of trl: event code
      
      
      % get lay
      load([info.savepath currsubj '_fieldtrip_layout.mat'])
      [dum laychanind] = match_str(nwaycomp.label,lay.label);
      
      % loop over comps
      for icomp = 1:ncomp
        figure('numbertitle','off','name',[currsubj '-' namesuff ' c' num2str(icomp)]);
        
        % gather
        A = comp{icomp}{1};
        B = comp{icomp}{2};
        if hashfa
          D = comp{icomp}{3};
        else
          C  = comp{icomp}{3};
          D = comp{icomp}{4};
        end
        if hasrpt
          if hashfa
            E = comp{icomp}{4};
          else
            E = comp{icomp}{5};
          end
        end
        if icomptrlload == 1
          Ec = E;
          E = abs(E);
        end
        
        % AMPCHAN/PHASCHAN
        for iset = 1:2
          subplot(3,6,(1:2) + (2*(iset-1)))
          hold on
          % plot outlines
          for ioutline = 1:numel(lay.outline)
            xline = lay.outline{ioutline}(:,1);
            yline = lay.outline{ioutline}(:,2);
            line(xline,yline,'linewidth',1,'linestyle','--')
          end
          % plot circles of A
          x = lay.pos(laychanind,1);
          y = lay.pos(laychanind,2);
          z = zeros(size(x));
          if iset == 1
            mag  = abs(A);
            phas = angle(A);
          else
            mag  = abs(B);
            phas = angle(B);
          end
          chansiz = lay.width(1)*1.5*mag;
          cfg = [];
          cfg.markercolor   = 'cparam';
          cfg.viewpoint     = [0 90];
          cfg.coordinates   = [x y z];
          cfg.cparam        = phas;
          cfg.renderbrain   = 'no';
          cfg.colormap      = 'hsv';
          cfg.clim          = [-pi pi];
          cfg.colorbar      = 'yes';
          cfg.axissqueeze   = 'yes';
          cfg.markersize    = num2cell(chansiz);
          roe_brainplot_chancircle(cfg);
          if iset == 1
            title('amp chan')
          else
            title('phas chan')
          end
        end
        
        % SPECTRAL PROFILES
        subplot(3,6,[5 6])
        plot(nwaycomp.phasfreq,D,'r')
        hold on
        if ~hashfa
          plot(nwaycomp.ampfreq,C,'b')
          legend({'phasfreq','ampfreq'})
        else
          legend('phasfreq')
        end
        axis tight
        xlabel('Hz')
        ylabel('loading')
        ylim = get(gca,'ylim');
        set(gca,'ylim',[min([0 ylim(1)]) ylim(2)])
        title('freq prof.')
        
        
        % TRIAL PROFILES
        subplot(3,6,7)
        plot(E);
        axis tight
        ylim = get(gca,'ylim');
        set(gca,'ylim',[min(0, min(ylim)) ylim(2)*1.05])
        xlabel('trial')
        ylabel('loading')
        % att/ign and mon/bin att
        subplot(3,6,[8 9])
        ind1 = att;
        ind2 = ign;
        m1 = mean(E(ind1));
        m2 = mean(E(ind2));
        sem1 = std(E(ind1)) ./ sqrt(sum(ind1));
        sem2 = std(E(ind2)) ./ sqrt(sum(ind2));
        ind3 = mon;
        ind4 = bin;
        m3 = mean(E(ind3));
        m4 = mean(E(ind4));
        sem3 = std(E(ind3)) ./ sqrt(sum(ind3));
        sem4 = std(E(ind4)) ./ sqrt(sum(ind4));
        errorbar([m1 m2; m3 m4]',[sem1 sem2; sem3 sem4]')
        ylim = get(gca,'ylim');
        set(gca,'ylim',[min(0,ylim(1)) ylim(2)*1.1]);
        set(gca,'xtick',[])
        set(gca,'ytick',[min(0,ylim(1)) ylim(2)*1.1])
        ylabel('loading')
        legend({['att-ign ' num2str(sum(ind1)) '/' num2str(sum(ind2))],['mon-bin ' num2str(sum(ind2)) '/' num2str(sum(ind3))]})
        [h p,ci,stats12] = ttest2(E(ind1),E(ind2));
        [h p,ci,stats34] = ttest2(E(ind3),E(ind4));
        title({['att/ign t = ' num2str(stats12.tstat,'%1.2f') '    ' 'mon/bin t = ' num2str(stats34.tstat,'%1.2f')]})
        % ipsi/contra lo/hi
        subplot(3,6,10)
        ind1 = ipsi;
        ind2 = contra;
        m1 = mean(E(ind1));
        m2 = mean(E(ind2));
        sem1 = std(E(ind1)) ./ sqrt(sum(ind1));
        sem2 = std(E(ind2)) ./ sqrt(sum(ind2));
        ind3 = lo;
        ind4 = hi;
        m3 = mean(E(ind3));
        m4 = mean(E(ind4));
        sem3 = std(E(ind3)) ./ sqrt(sum(ind3));
        sem4 = std(E(ind4)) ./ sqrt(sum(ind4));
        errorbar([m1 m2; m3 m4]',[sem1 sem2; sem3 sem4]')
        ylim = get(gca,'ylim');
        set(gca,'ylim',[min(0,ylim(1)) ylim(2)*1.1]);
        set(gca,'xtick',[])
        set(gca,'ytick',[min(0,ylim(1)) ylim(2)*1.1])
        ylabel('loading')
        legend({['ipsi-contra ' num2str(sum(ind1)) '/' num2str(sum(ind2))],['lo-hi ' num2str(sum(ind2)) '/' num2str(sum(ind3))]})
        [h p,ci,stats12] = ttest2(E(ind1),E(ind2));
        [h p,ci,stats34] = ttest2(E(ind3),E(ind4));
        title({['ipsi/contra t = ' num2str(stats12.tstat,'%1.2f') '    ' 'lo/hi t = ' num2str(stats34.tstat,'%1.2f')]})
        %
        subplot(3,6,[11 12])
        ind1 = att & ipsi;
        ind2 = ign & ipsi;
        m1 = mean(E(ind1));
        m2 = mean(E(ind2));
        sem1 = std(E(ind1)) ./ sqrt(sum(ind1));
        sem2 = std(E(ind2)) ./ sqrt(sum(ind2));
        ind3 = att & contra;
        ind4 = ign & contra;
        m3 = mean(E(ind3));
        m4 = mean(E(ind4));
        sem3 = std(E(ind3)) ./ sqrt(sum(ind3));
        sem4 = std(E(ind4)) ./ sqrt(sum(ind4));
        errorbar([m1 m2; m3 m4]',[sem1 sem2; sem3 sem4]')
        ylim = get(gca,'ylim');
        set(gca,'ylim',[min(0,ylim(1)) ylim(2)*1.1]);
        set(gca,'xtick',[])
        set(gca,'ytick',[min(0,ylim(1)) ylim(2)*1.1])
        ylabel('loading')
        legend({['att-ign ipsi ' num2str(sum(ind1)) '/' num2str(sum(ind2))],['att-ign contra ' num2str(sum(ind2)) '/' num2str(sum(ind3))]})
        [h p,ci,stats12] = ttest2(E(ind1),E(ind2));
        [h p,ci,stats34] = ttest2(E(ind3),E(ind4));
        title({['att/ign ipsi t = ' num2str(stats12.tstat,'%1.2f') '    ' 'att/ign contra t = ' num2str(stats34.tstat,'%1.2f')]})
        
        
        % att/ign and mon/bin att
        subplot(3,6,13)
        ind1 = att;
        ind2 = ign;
        compass(Ec(ind1))
        hold on
        compass(Ec(ind2),'r')
        [p table12]      = circ_wwtest(exp(1i.*angle(Ec(ind1))),exp(1i.*angle(Ec(ind2))));
        title({['att/ign f = ' num2str(table12{2,5},'%1.2f')]})
        subplot(3,6,14)
        ind3 = mon;
        ind4 = bin;
        compass(Ec(ind3))
        hold on
        compass(Ec(ind4),'r')
        [p table34]      = circ_wwtest(exp(1i.*angle(Ec(ind3))),exp(1i.*angle(Ec(ind4))));
        title({['mon/bin f = ' num2str(table34{2,5},'%1.2f')]})
        % ipsi/contra lo/hi
        subplot(3,6,15)
        ind1 = ipsi;
        ind2 = contra;
        compass(Ec(ind1))
        hold on
        compass(Ec(ind2),'r')
        [p table12]      = circ_wwtest(exp(1i.*angle(Ec(ind1))),exp(1i.*angle(Ec(ind2))));
        title({['ipsi/contra f = ' num2str(table12{2,5},'%1.2f')]})
        subplot(3,6,16)
        ind3 = lo;
        ind4 = hi;
        compass(Ec(ind3))
        hold on
        compass(Ec(ind4),'r')
        [p table34]      = circ_wwtest(exp(1i.*angle(Ec(ind3))),exp(1i.*angle(Ec(ind4))));
        title({['lo/hi f = ' num2str(table34{2,5},'%1.2f')]})
        %
        subplot(3,6,17)
        ind1 = att & ipsi;
        ind2 = ign & ipsi;
        compass(Ec(ind1))
        hold on
        compass(Ec(ind2),'r')
        [p table12]      = circ_wwtest(exp(1i.*angle(Ec(ind1))),exp(1i.*angle(Ec(ind2))));
        title({['att/ign ipsi f = ' num2str(table12{2,5},'%1.2f')]})
        subplot(3,6,18)
        ind3 = att & contra;
        ind4 = ign & contra;
        compass(Ec(ind3))
        hold on
        compass(Ec(ind4),'r')
        [p table34]      = circ_wwtest(exp(1i.*angle(Ec(ind3))),exp(1i.*angle(Ec(ind4))));
        title({['att/ign contra f = ' num2str(table34{2,5},'%1.2f')]})
        
             
             
      end % icomp
      
      
      
    end
  end
end



























