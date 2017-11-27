function rmr_predfaceval_getspace




% fetch info
info = rmr_predfaceval_info;

%info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

% set name suffix
%fnnamesuffix = '4to40hz_3cyc_CTI';
fnnamesuffix = '6to40hz_3cyc_postface';

% nway settings
nwayalg     = 'spacefsp'; % 'spacefsp'
nwaynmethod = 'ncomp10'; %'ncomp15sh'; % 'splitrel' 'ncomp10sh'
nwaynrand    = 20;
nwayconvcrit = 1e-8;
normmethod   = 'coh'; % 'coh'   'avgoverepoch' '16throotpower'
nwaysplit    = 'oddeventrials'; %  oddevenspikes  oddeventrials
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];

% read data and get fourier, 1 session
distcompsys = 'torque';
for      isubj = 2:3%     :numel(info.subj)
  
  dosplitrel = strcmp(nwaynmethod,'splitrel') || strcmp(nwaynmethod(end-1:end),'sr');
  if dosplitrel
    nwayadd = [nwayadd '_' 'split' nwaysplit];
  end
  
  % set currs
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj ', ' fnnamesuffix ', ' nwayadd])
  
  % set nway filename and continue if it doesn't exist yet
  nwayfn = [info.savepath currsubj '_' fnnamesuffix '_' nwayadd '.mat'];
  if ~exist(nwayfn,'file')
    
    
    % set fns and load data
    fourfn = [info.savepath currsubj '_' 'fourier' '_' fnnamesuffix '.mat'];
    datafn = [info.savepath currsubj '_' 'dataetc' '_' fnnamesuffix '.mat'];
    load(datafn)
    data = datasel;
    
    % create normed fourier
    fourndfn = [info.scratchpath currsubj '_' 'fourier' '_' fnnamesuffix '_' normmethod '.mat'];
    load(fourfn)
    if strncmp(fliplr(normmethod),fliplr('throotpower'),11)
      [fourier,scaling] = roe_fouriernormalize(fourier,normmethod);
    else
      fourier = roe_fouriernormalize(fourier,normmethod);
    end
    save(fourndfn,'fourier','-v7.3')
    clear fourier
    
    % create split and save
    if dosplitrel
      switch nwaysplit
        
        case 'oddeventrials'
          % load from disk
          load(fourndfn)
          
          % set split indices
          ntrial = size(fourier,3);
          trlindprt1 = 1:2:ntrial;
          trlindprt2 = 2:2:ntrial;
          
          % save splits
          fourierfull = fourier;
          % set fns
          fourprt1fn = [info.scratchpath currsubj '_' 'fourier' '_' fnnamesuffix '_' normmethod '_' nwaysplit '_' 'prt1' '.mat'];
          fourprt2fn = [info.scratchpath currsubj '_' 'fourier' '_' fnnamesuffix '_' normmethod '_' nwaysplit '_' 'prt2' '.mat'];
          % save split1
          fourier = fourierfull(:,:,trlindprt1,:);
          save(fourprt1fn,'fourier','-v7.3')
          clear fourier
          % save split2
          fourier = fourierfull(:,:,trlindprt2,:);
          save(fourprt2fn,'fourier','-v7.3')
          clear fourier
          clear fourierfull
          % save fns for fourierdata
          fourierpart = {fourprt1fn,fourprt2fn};
          
        otherwise
          error('nwaysplit not supported')
      end
    end
    
    
    % build fourierdata
    fourierdata = [];
    fourierdata.fourier        = fourndfn;
    if dosplitrel
      fourierdata.fourierpart  = fourierpart;
    end
    fourierdata.freq           = freqoi;
    fourierdata.label          = data.label;
    fourierdata.dimord         = 'chan_freq_epoch_tap';
    fourierdata.trialinfo      = data.trialinfo;
    fourierdata.cfg            = data.cfg;
    clear fourier
    
    % nway settings
    cfg = [];
    cfg.model              = nwayalg;
    cfg.datparam           = 'fourier';
    cfg.Dmode              = 'identity';
    cfg.numiter            = 3000;
    cfg.convcrit           = nwayconvcrit;
    cfg.randstart          = nwaynrand;
    cfg.outputfile         = nwayfn;
    cfg.checkpointpath     = [info.scratchpath 'CHECKPOINT' '_' 'predfaceval' '_' nwaynmethod '_' normmethod '_'];
    % qsub options
    switch distcompsys
      case 'torque'
        cfg.distcomp.system          = 'torque';
        cfg.distcomp.timreq          = 60*60*9; %
        cfg.distcomp.matlabcmd       = '/opt/matlab/2015a/bin/matlab';
        cfg.distcomp.torquequeue     = 'hotel';
        cfg.distcomp.torquestack     = [];
        cfg.distcomp.inputsaveprefix = '/oasis/tscc/scratch/roevdmei/nwaytemp_';
        cfg.distcomp.qsuboptions     = ' -k oe ';
      case 'matlabpct'
        cfg.distcomp.system          = 'matlabpct';
        if info.ontscc
          cfg.distcomp.mpctpoolsize  = str2double(getenv('PBS_NUM_PPN'));
          cfg.distcomp.mpctcluster   = parcluster('local');
          scratchdir                 = [info.scratchpath 'matlabpool' num2str(round(sum(clock*1e6))) '/'];
          mkdir(scratchdir);
          cfg.distcomp.mpctcluster.JobStorageLocation = scratchdir;
        end
    end
    % ncomp options
    if strcmp(nwaynmethod(1:5),'ncomp') && ~strcmp(nwaynmethod(end-1:end),'sr')
      cfg.ncomp = str2double(nwaynmethod(6:end));
    elseif strcmp(nwaynmethod(1:5),'ncomp') && strcmp(nwaynmethod(end-1:end),'sr')
      cfg.ncompestrandstart   = nwaynrand;
      cfg.ncompestsrdatparam  = 'fourierpart';
      cfg.ncompest            = 'splitrel';
      cfg.ncompeststart       = str2double(nwaynmethod(6:end-2));
      cfg.ncompeststep        = 1;
      cfg.ncompestend         = str2double(nwaynmethod(6:end-2));
      cfg.ncompestsrcritval   = [NaN NaN 0 NaN 0];
      cfg.ncompestsrcritjudge = 'meanoversplitscons';
    elseif strcmp(nwaynmethod,'splitrel')
      cfg.ncompestrandstart   = nwaynrand;
      cfg.ncompestsrdatparam  = 'fourierpart';
      cfg.ncompest            = 'splitrel';
      cfg.ncompeststart       = 10;
      cfg.ncompeststep        = 10;
      cfg.ncompestend         = 80;
      cfg.ncompestsrcritval   = [.7 .7 0 .7 0];
      cfg.ncompestsrcritjudge = 'meanoversplitscons';
    else
      error('...')
    end
    % go for it
    nwaycomp = nd_nwaydecomposition(cfg,fourierdata);
    
    % save on oasis just in case
    if info.ontscc
      save([info.scratchpath currsess '_' fnnamesuffix '_' nwayadd '.mat'],'nwaycomp','-v7.3');
      % remove scratch if ontscc
      if strcmp(distcompsys,'matlabpct') && exist(cfg.distcomp.mpctcluster.JobStorageLocation,'dir')
        delete([cfg.distcomp.mpctcluster.JobStorageLocation '/' 'matlab_metadata.mat']) % should be the only file present
        rmdir(cfg.distcomp.mpctcluster.JobStorageLocation) % only works if folder is empty, which it now should be
      end
    end
    
    % clear remaining stuff
    clear fourier nwaycomp fourierdata
    
  end % ~exist(nwayfn,'file')
end % isubj





% get out of here
exit












function playground






% fetch info
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';
info.savepath = '/Users/roemer/Work/Data/tmpdata/predfaceval/';

% set name suffix
fnnamesuffix = '4to40hz_3cyc_CTI';
%fnnamesuffix = '6to40hz_3cyc_postface';

% nway settings
nwayalg     = 'spacefsp'; % 'spacefsp'
nwaynmethod = 'ncomp10'; %'ncomp15sh'; % 'splitrel' 'ncomp10sh'
nwaynrand    = 20;
nwayconvcrit = 1e-8;
normmethod   = 'none'; % 'coh'   'avgoverepoch' '16throotpower'
nwaysplit    = 'oddeventrials'; %  oddevenspikes  oddeventrials
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];

% read data and get fourier, 1 session
distcompsys = 'torque';
for      isubj = 2:3%     :numel(info.subj)
  
  dosplitrel = strcmp(nwaynmethod,'splitrel') || strcmp(nwaynmethod(end-1:end),'sr');
  if dosplitrel
    nwayadd = [nwayadd '_' 'split' nwaysplit];
  end
  
  % set currs
  currsubj = info.subj{isubj};
  
  % set fns and load data
  fourfn = [info.savepath currsubj '_' 'fourier' '_' fnnamesuffix '.mat'];
  datafn = [info.savepath currsubj '_' 'dataetc' '_' fnnamesuffix '.mat'];
  load(datafn)
  data = datasel;
  
  % set nway filename and continue if it doesn't exist yet
  nwayfn = [info.savepath currsubj '_' fnnamesuffix '_' nwayadd '.mat'];
  if exist(nwayfn,'file')
    load(nwayfn)
    
    % plot decomp stats
    cfg = [];
    cfg.savefig    = 'no';
    cfg.figprefix  = 'paco';
    cfg.figvisible = 'on';
    roe_nwaystatplot(cfg,nwaycomp) % OLD SHITTY FUNCTION, TO INVESTIGATE RANDOM STARTS AND SPLIT RELIABILITY SEE TUTORIALS
    
    % plot congrcumul
    figure('numbertitle','off','name','congrcumul')
    congrcumul = nwaycomp.randomstat.congrcumul;
    nrand = size(congrcumul,1);
    plottitle = {'A','B','C','S'};
    for iparam = 1:4
      subplot(2,2,iparam)
      plot(2:nrand,squeeze(congrcumul(2:end,:,iparam)))
      xlabel('1:N random starts')
      ylabel('congruence')
      set(gca,'xlim',[2 nrand],'ylim',[0.5 1])
      title(plottitle{iparam})
    end
    
    
    % set basics
    comp   = nwaycomp.comp;
    ncomp  = numel(comp);
    nchan  = size(comp{1}{1},1);
    nfreq  = size(comp{1}{2},1);
    ntrial = size(comp{1}{3},1);
    
    
    % gather A,B,C,S
    A = zeros(nchan,ncomp);
    B = zeros(nfreq,ncomp);
    C = zeros(ntrial,ncomp);
    %S = zeros(nchan,ncomp);
    L = zeros(nchan,ncomp);
    for icomp = 1:ncomp
      A(:,icomp) = comp{icomp}{1};
      B(:,icomp) = comp{icomp}{2};
      C(:,icomp) = comp{icomp}{3};
      %S(:,icomp) = comp{icomp}{4};
      L(:,icomp) = comp{icomp}{4}(:,max(comp{icomp}{2})==comp{icomp}{2});
    end
    
    % get layout
    for itrial = 1:numel(data.time);
      data.trial{itrial} = rand(numel(data.label),numel(data.time{itrial}));
    end
    cfg = [];
    cfg.layout = 'ordered';
    lay = ft_prepare_layout(cfg,data);
    ind = [];
    type = {'RAM','LAM','AM'};
    for itype = 1:numel(type)
      ind = [ind; find(strncmp(lay.label,type{itype},numel(type{itype})))];
    end
    ind = [ind; setdiff(1:nchan+2,ind)'];
    lay.label  = lay.label(ind);
    
    % sel
    [dum, laychanind] = match_str(nwaycomp.label,lay.label);

    
    % per comp
    for icomp = 1:ncomp
      figure('numbertitle','off','name',['comp' num2str(icomp)]);
      % plot A
      subplot(2,3,1)
      % plot circles of A
      x = lay.pos(laychanind,1);
      y = lay.pos(laychanind,2);
      z = zeros(size(x));
      mag  = abs(A(:,icomp));
      phas = angle(exp(1i*L(:,icomp) .* 2*pi));
      chansiz = lay.width(1)*500*mag;
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
      cfg.labelsize     = 6;
      roe_brainplot_chancircle(cfg);
      title('A and L at peak freq')
      
      % plot B
      subplot(2,3,2)
      plot(nwaycomp.freq,B(:,icomp));
      axis tight
      ylim = get(gca,'ylim');
      set(gca,'ylim',[0 ylim(2)*1.05])
      xlabel('frequency (index)')
      ylabel('loading')
      title('B')
      
%       % plot S
%       subplot(2,3,3)
%       hold on
%       [dum sortind] = sort(A(:,icomp),'descend');
%       plot(S(sortind(1:10),icomp));
%       plot(diff(S(sortind(1:10),icomp)),'r');
%       ylim = get(gca,'ylim');
%       axis tight
%       set(gca,'ylim',[-.025 .025])
%       xlabel('sorted units')
%       ylabel('loading')
%       title('S')
      
      % plot C
      subplot(2,3,4)
      [dum sortind] = sort(C(:,icomp),'descend');
      [ax h1 h1] = plotyy(1:ntrial,C(	sortind,icomp),1:ntrial,nwaycomp.trialinfo(sortind,[1 2]));
      ylim = get(ax(1),'ylim');
      set(ax(1),'ylim',[0 ylim(2)*1.05],'xlim',[1 ntrial])
      ylim = get(ax(2),'ylim');
      set(ax(2),'ylim',[1 ylim(2)*1.05],'xlim',[1 ntrial])
      xlabel('trial')
      ylabel('loading')
      
      % C t-test: pred vs unpred, fear vs neutral
      subplot(2,3,5)
      if strncmp(fliplr(fnnamesuffix),'ITC',3)
        ind1 = nwaycomp.trialinfo(:,1)==1; % pred vs unpred
        ind2 = ~ind1;
        m1 = mean(C(ind1,icomp));
        m2 = mean(C(ind2,icomp));
        sem1 = std(C(ind1,icomp)) ./ sqrt(sum(ind1));
        sem2 = std(C(ind2,icomp)) ./ sqrt(sum(ind2));
        errorbar([m1 m2],[sem1 sem2])
        ylim = get(gca,'ylim');
        set(gca,'ylim',[min(0,ylim(1)) ylim(2)*1.1]);
        set(gca,'xtick',[])
        set(gca,'ytick',[min(0,ylim(1)) ylim(2)*1.1])
        ylabel('loading')
        legend({['pred-unpred ' num2str(sum(ind1)) '/' num2str(sum(ind2))]})
        [h p,ci,stats12] = ttest2(C(ind1,icomp),C(ind2,icomp));
        title(['pred/unpred t = ' num2str(stats12.tstat,'%1.2f')])
      else
        ind1 = nwaycomp.trialinfo(:,1)==1; % pred vs unpred
        ind2 = ~ind1;
        m1 = mean(C(ind1,icomp));
        m2 = mean(C(ind2,icomp));
        sem1 = std(C(ind1,icomp)) ./ sqrt(sum(ind1));
        sem2 = std(C(ind2,icomp)) ./ sqrt(sum(ind2));
        ind3 = nwaycomp.trialinfo(:,2)==1; % fear vs neut
        ind4 = ~ind3;
        m3 = mean(C(ind3,icomp));
        m4 = mean(C(ind4,icomp));
        sem3 = std(C(ind3,icomp)) ./ sqrt(sum(ind3));
        sem4 = std(C(ind4,icomp)) ./ sqrt(sum(ind4));
        errorbar([m1 m2; m3 m4]',[sem1 sem2; sem3 sem4]')
        ylim = get(gca,'ylim');
        set(gca,'ylim',[min(0,ylim(1)) ylim(2)*1.1]);
        set(gca,'xtick',[])
        set(gca,'ytick',[min(0,ylim(1)) ylim(2)*1.1])
        ylabel('loading')
        legend({['pred-unpred ' num2str(sum(ind1)) '/' num2str(sum(ind2))],['fear-neut ' num2str(sum(ind2)) '/' num2str(sum(ind3))]})
        [h p,ci,stats12] = ttest2(C(ind1,icomp),C(ind2,icomp));
        [h p,ci,stats34] = ttest2(C(ind3,icomp),C(ind4,icomp));
        % anova
        %[p,tabel,stats] = anovan(E,{pred,fear});
        title(['pred/unpred t = ' num2str(stats12.tstat,'%1.2f') '    ' 'fear/neut t = ' num2str(stats34.tstat,'%1.2f')])
      end
      
%       if strncmp(fliplr(fnnamesuffix),'ITC',3)
%         ind1 = nwaycomp.trialinfo(:,1)==1; % pred vs unpred
%       else
%         ind1 = nwaycomp.trialinfo(:,2)==1; % fear vs neut
%       end
%       ind2 = ~ind1;
%       m1 = mean(C(ind1,icomp));
%       m2 = mean(C(ind2,icomp));
%       sem1 = std(C(ind1,icomp)) ./ sqrt(sum(ind1));
%       sem2 = std(C(ind2,icomp)) ./ sqrt(sum(ind2));
%       errorbar([m1 m2],[sem1 sem2])
%       ylim = get(gca,'ylim');
%       set(gca,'ylim',[0 ylim(2)*1.1]);
%       set(gca,'xtick',[1 2])
%       set(gca,'ytick',[0 ylim(2)*1.1])
%       if strncmp(fliplr(fnnamesuffix),'ITC',3)
%         set(gca,'xticklabel',{'pred','unpred'})
%       else
%         set(gca,'xticklabel',{'fearful','neutral'})
%       end
%       ylabel('loading')
%       [h p,ci,stats] = ttest2(C(ind1,icomp),C(ind2,icomp));
%       title(['t = ' num2str(stats.tstat,'%1.2f')])
      
      % C
      %subplot(2,3,6)
    enc
    
    
  end % exist
  
  
end







































