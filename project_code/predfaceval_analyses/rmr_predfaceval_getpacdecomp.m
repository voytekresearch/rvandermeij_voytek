function rmr_predfaceval_getpacdecomp





% fetch info
info = rmr_predfaceval_info;

%info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

% design
centerampfreq = 'yes';

% nway
nwayalg     = 'parafac';
nwaynmethod = 'corcondiag'; % 'corcondiag'  'splithalf'
nwayadd     = 'rnd100_conv1e-6'; %
%
%comptrlload = 0;

for     isubj = 1:numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
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
  
  
  % for artifacts been detected
  if ~info.(currsubj).trlartfctflg
    disp(['artifacts not yet detected for ' currsubj])
    continue
  end
  
  %
  datfn = [];
  datfn{1} = [info.savepath currsubj '_' 'wplf_4to60'           '_'         'CTI' '_' 'centerae' centerampfreq  '.mat'];
  datfn{2} = [info.savepath currsubj '_' 'wplf_4to60to70to130'  '_'         'CTI' '_' 'centerae' centerampfreq  '.mat'];
  datfn{3} = [info.savepath currsubj '_' 'wplf_4to60toHFA'      '_'         'CTI' '_' 'centerae' centerampfreq  '.mat'];
  datfn{4} = [info.savepath currsubj '_' 'wplf_6to60'           '_' 'postface500' '_' 'centerae' centerampfreq  '.mat'];
  datfn{5} = [info.savepath currsubj '_' 'wplf_6to60to70to130'  '_' 'postface500' '_' 'centerae' centerampfreq  '.mat'];
  datfn{6} = [info.savepath currsubj '_' 'wplf_6to60toHFA'      '_' 'postface500' '_' 'centerae' centerampfreq  '.mat'];
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
      if strcmp(nwaynmethod,'corcondiag')
        load(datfn{idataset})
        wplf = squeeze(wplfdata.wplf);
        if false%ndims(wplf) == 4
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
        if false%ndims(wplf) == 4
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
        cfg.ncompestshcritval  = [.7 .7 .7 .7];     % splithalf criterion for each of the parameters
      else
        cfg.complexdims        = [1 1 0];
        cfg.ncompestshcritval  = [.7 .7 .7];     % splithalf criterion for each of the parameters
      end
      if false;%ndims(wplf) == 4
        % qsub options
        cfg.distcomp.system          = 'torque';
        cfg.distcomp.timreq          = 60*60*24*.5; % 0.5 days
        %cfg.distcomp.inputsaveprefix = '/home/roevdmei/qsubtemp/nwaytemp_';
        cfg.distcomp.matlabcmd       = '/opt/matlab/2014b/bin/matlab';
        cfg.distcomp.torquequeue     = 'hotel';
      end
      %
      nwaycomp = nd_nwaydecomposition(cfg,wplfdata);
    end
    
    
    
    
    
    %%% OBTAIN TRIAL LOADINGS
    for comptrlload = 1
      if ~exist([nwayfn{idataset}(1:end-4) '_' 'trialloadings' '_' 'comp' num2str(comptrlload) '.mat'],'file')
        load(datfn{idataset})
        load(nwayfn{idataset})
        
        % fetch and preprocess data
        % combine data over sessions
        currfndat = [info.savepath currsubj '_' 'allsess_bs540lp_detrend_0mspre_600mspost' '.mat'];
        if ~exist(currfndat,'file')
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
          save(currfndat,'data','-v7.3');
        else
          load(currfndat)
        end
        
        % cut out data
        cfg = [];
        cfg.toilim = wplfdata.cfg.previous.previous.previous.toilim;
        datasel = ft_redefinetrial(cfg,data);
        cfg = [];
        cfg.demean  = 'yes';
        cfg.detrend = 'yes';
        cfg.trials  = 'all';
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
info = rmr_predfaceval_info;
info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';


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
  for iperiod = 1:2
    datfn = [];
    if iperiod == 1
      datfn{1} = [info.savepath currsubj '_' 'wplf_4to60toHFA' '_'           'CTI' '_' 'centerae' centerampfreq '.mat'];
      %datfn{2} = [info.savepath currsubj '_' 'wplf_4to150'     '_'           'CTI' '_' 'centerae' centerampfreq '.mat'];
    elseif iperiod == 2
      datfn{1} = [info.savepath currsubj '_' 'wplf_6to60toHFA' '_'   'postface500' '_' 'centerae' centerampfreq '.mat'];
      %datfn{2} = [info.savepath currsubj '_' 'wplf_6to150'     '_'   'postface500' '_' 'centerae' centerampfreq '.mat'];
    end
    %
    for idataset = 1:numel(datfn)
      for icomptrlload = 1
        
        nwayfn = [datfn{idataset}(1:end-4) '_' nwayalg '_' nwaynmethod '_' nwayadd '_' 'trialloadings'  '_' 'comp' num2str(icomptrlload) '.mat'];
        if ~exist(nwayfn,'file')
          continue
        end
        
        load(nwayfn)
        switch idataset
          case {1}
            namesuff = 'hfa';
          case {2}
            namesuff = 'osc';
        end
        switch iperiod
          case {1}
            namesuff = [namesuff '_' 'cti'];
          case {2}
            namesuff = [namesuff '_' 'pf'];
        end
        
        
        %               % plot decomp stats
        %               cfg = [];
        %               cfg.savefig    = 'no';
        %               cfg.figvisible = 'on';
        %               roe_nwaystatplot(cfg,nwaycomp)
        
        
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
        
        % for CTI: pred/unpred
        % for postface: fear/neut
        % for postface: pred/unpred
        % for postface: fear/neut and pred/unpred
        % get them all, and select later
        valencetc = nwaycomp.trialinfo(:,2); % fearful (1) or neutral (2) face trial
        predtc    = nwaycomp.trialinfo(:,1); % pred (1) or unpred (2) trial
        pred   = predtc == 1;
        unpred = predtc == 2;
        fear   = valencetc == 1;
        neut   = valencetc == 2;
        
        % get layout
        cfg = [];
        cfg.layout = 'ordered';
        lay = ft_prepare_layout(cfg,nwaycomp);       % reorder channels in lay
        ind = [];
        type = {'RAM','LAM','AM'};
        for itype = 1:numel(type)
          ind = [ind; find(strncmp(lay.label,type{itype},numel(type{itype})))];
        end
        ind = [ind; setdiff(1:nchan+2,ind)'];
        lay.label  = lay.label(ind);
        
        % sel
        [dum, laychanind] = match_str(nwaycomp.label,lay.label);
        
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
            E = abs(E);
          end
          
          % AMPCHAN/PHASCHAN
          for iset = 1:2
            subplot(2,3,iset)
            hold on
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
            cfg.label         = nwaycomp.label;
            cfg.labelcolor    = 'k';
            cfg.labelsize     = 8;
            roe_brainplot_chancircle(cfg);
            if iset == 1
              title('amp chan')
            else
              title('phas chan')
            end
          end
          
          % SPECTRAL PROFILES
          subplot(2,3,3)
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
          subplot(2,3,4)
          plot(E);
          axis tight
          ylim = get(gca,'ylim');
          set(gca,'ylim',[min(0, min(ylim)) ylim(2)*1.05])
          xlabel('trial')
          ylabel('loading')
          %
          if iperiod == 1
            subplot(2,3,5)
            ind1 = pred;
            ind2 = unpred;
            m1 = mean(E(ind1));
            m2 = mean(E(ind2));
            sem1 = std(E(ind1)) ./ sqrt(sum(ind1));
            sem2 = std(E(ind2)) ./ sqrt(sum(ind2));
            errorbar([m1 m2],[sem1 sem2])
            ylim = get(gca,'ylim');
            set(gca,'ylim',[min(0,ylim(1)) ylim(2)*1.1]);
            set(gca,'xtick',[])
            set(gca,'ytick',[min(0,ylim(1)) ylim(2)*1.1])
            ylabel('loading')
            legend({['pred-unpred ' num2str(sum(ind1)) '/' num2str(sum(ind2))]})
            [h p,ci,stats12] = ttest2(E(ind1),E(ind2));
            title(['pred/unpred t = ' num2str(stats12.tstat,'%1.2f')])
          elseif iperiod == 2
            subplot(2,3,5)
            ind1 = pred;
            ind2 = unpred;
            m1 = mean(E(ind1));
            m2 = mean(E(ind2));
            sem1 = std(E(ind1)) ./ sqrt(sum(ind1));
            sem2 = std(E(ind2)) ./ sqrt(sum(ind2));
            ind3 = fear;
            ind4 = neut;
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
            legend({['pred-unpred ' num2str(sum(ind1)) '/' num2str(sum(ind2))],['fear-neut ' num2str(sum(ind2)) '/' num2str(sum(ind3))]})
            [h p,ci,stats12] = ttest2(E(ind1),E(ind2));
            [h p,ci,stats34] = ttest2(E(ind3),E(ind4));
            % anova
            %[p,tabel,stats] = anovan(E,{pred,fear});
            title(['pred/unpred t = ' num2str(stats12.tstat,'%1.2f') '    ' 'fear/neut t = ' num2str(stats34.tstat,'%1.2f')])
          end
          
          
        end % icomp
        
        
        
      end % comptrlload
    end % idataset
  end % iperiod
end















