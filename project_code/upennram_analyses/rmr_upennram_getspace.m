function rmr_upennram_getspace



% get info
info = rmr_upennram_info;
info.subj = info.subjselmains;


% set name suffix
fnnamesuffix = '3to40hz_3cyc';

% nway settings
nwayalg      = 'spacetime'; % 'spacefsp'
nwaynmethod  = 'ncomp20sr'; %'ncomp15sr'; % 'splitrel'
nwaynrand    = 50; % 'rnd40_conv1e-6'
nwayconvcrit = 1e-8;
normmethod   = 'none'; % coh none 16throotpower
nwaysplit    = 'oddeventrials'; %  oddeventrials
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];

% read data and get fourier
distcompsys = 'torque';
dosplitrel = strcmp(nwaynmethod,'splitrel') || strcmp(nwaynmethod(end-1:end),'sr');
if dosplitrel
  nwayadd = [nwayadd '_' 'split' nwaysplit];
end
for     isubj = 1 % :numel(info.subj)
  
  % set
  currsubj   = info.subj{isubj};
  disp(['working on ' currsubj ', ' fnnamesuffix ', ' nwayadd])
  
  % set nway filename and continue if it doesn't exist yet
  nwayfn = [info.savepath currsubj '_' fnnamesuffix '_' nwayadd '.mat'];
  if ~exist(nwayfn,'file')
    
    % set fns and load data
    fourfn = [info.savepath currsubj '_' 'fourier' '_' fnnamesuffix '.mat'];
    datafn = [info.savepath currsubj '_' 'dataetc' '_' fnnamesuffix '.mat'];
    load(datafn)
    
    % create normed fourier
    fourndfn = [info.scratchpath currsubj '_' 'fourier' '_' fnnamesuffix '_' normmethod '.mat'];
    load(fourfn)
    fourier = rmr_fouriernormalize(fourier,normmethod);
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
    cfg.checkpointpath     = ['/oasis/tscc/scratch/roevdmei/' 'CHECKPOINT' '_' 'upennram' '_' currsubj '_' nwaynmethod '_' normmethod '_'];
    % qsub options
    switch distcompsys
      case 'torque'
        cfg.distcomp.system          = 'torque';
        cfg.distcomp.timreq          = 60*60*30; %
        cfg.distcomp.matlabcmd       = '/opt/matlab/2015a/bin/matlab';
        cfg.distcomp.torquequeue     = 'hotel';
        cfg.distcomp.torquestack     = 1;
        %cfg.distcomp.inputsaveprefix = '/oasis/tscc/scratch/roevdmei/nwaytemp_';
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
      cfg.ncompeststart       = 5;
      cfg.ncompeststep        = 5;
      cfg.ncompestend         = 50;
      cfg.ncompestsrcritval   = [.7 .7 0 .7 0];
      cfg.ncompestsrcritjudge = 'meanoversplitscons';
    else
      error('...')
    end
    % go for it
    nwaycomp = nd_nwaydecomposition(cfg,fourierdata);
    
    % save on oasis just in case
    if info.ontscc
      save([info.scratchpath currsubj '_' fnnamesuffix '_' nwayadd '.mat'],'nwaycomp','-v7.3');
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



% get info
info = rmr_upennram_info;
info.subj = info.subjselmains;
%info.savepath = '/Users/roemervandermeij/Work/Data/tmpdata/upennram/';


% set name suffix
fnnamesuffix = '3to40hz_3cyc';

% nway settings
nwayalg      = 'spacetime'; % 'spacefsp'
nwaynmethod  = 'ncomp20sr'; %'ncomp15sr'; % 'splitrel'
nwaynrand    = 50; % 'rnd40_conv1e-6'
nwayconvcrit = 1e-8;
normmethod   = 'none'; % coh none 16throotpower
nwaysplit    = 'oddeventrials'; %  oddeventrials
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];

% read data and get fourier
if (strcmp(nwaynmethod,'splitrel') || strcmp(nwaynmethod(end-1:end),'sr'))
  nwayadd = [nwayadd '_' 'split' nwaysplit];
end
for     isubj = 1 % :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  
  
  % set nway filename and continue if it doesn't exist yet
  nwayfn = [info.savepath currsubj '_' fnnamesuffix '_' nwayadd '.mat'];
  if exist(nwayfn,'file')
    load(nwayfn)
    
    % plot decomp stats
    cfg = [];
    cfg.savefig    = 'no';
    cfg.figvisible = 'on';
    roe_nwaystatplot(cfg,nwaycomp)
    
    % plot congrcumul
    figure('numbertitle','off','name','congrcumul')
    congrcumul = nwaycomp.randomstat.congrcumul;
    %congrcumul = nwaycomp.splitrelstat.allrandomstatfull{isplit}.congrcumul;
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
    
    % gather A,B,C,S,L
    A = zeros(nchan,ncomp);
    B = zeros(nfreq,ncomp);
    C = zeros(ntrial,ncomp);
    switch nwaycomp.cfg.model
      case 'spacefsp'
        L = zeros(nchan,ncomp);
      case 'spacetime'
        S = zeros(nchan,ncomp);
    end
    for icomp = 1:ncomp
      A(:,icomp) = comp{icomp}{1};
      B(:,icomp) = comp{icomp}{2};
      C(:,icomp) = comp{icomp}{3};
      switch nwaycomp.cfg.model
        case 'spacefsp'
          L(:,icomp) = comp{icomp}{4}(:,max(comp{icomp}{2})==comp{icomp}{2});
        case 'spacetime'
          S(:,icomp) = comp{icomp}{4};
      end
    end
    
    % fetch eleccoords
    coords = rmr_upennram_fetcheleccoords(currsubj,nwaycomp.label);
    
    % per comp
    for icomp = 1:ncomp
      figure('numbertitle','off','name',['comp' num2str(icomp)]);
      % plot A
      subplot(3,4,[1 2 3 5 6 7 9 10 11])
      % plot circles of A
      mag     = abs(A(:,icomp));
      chansiz = 100*mag;
      switch nwaycomp.cfg.model
        case 'spacefsp'
          colparam  = angle(exp(1i*L(:,icomp) .* 2*pi));
          plotclim  = [-pi pi];
          plottitle = 'A and L at peak freq';
          colmap    = 'hsv';
        case 'spacetime'
          [dum maxind] = max(mag);
          colparam  = S(:,icomp)-S(maxind,icomp);
          plotclim  = [-0.05 0.05];
          plottitle = 'A with S';
          colmap    = 'parula';
      end
      cfg = [];
      cfg.markercolor   = 'cparam';
      cfg.viewpoint     = [90 0];
      cfg.coordinates   = coords;
      cfg.cparam        = colparam;
      cfg.renderbrain   = 'yes';
      cfg.template      = 'surface_white_both_reduced0.05.mat';
      cfg.alpha         = 0.1;
      cfg.colormap      = colmap;
      cfg.clim          = plotclim;
      cfg.colorbar      = 'yes';
      cfg.axissqueeze   = 'yes';
      cfg.markersize    = num2cell(chansiz);
      cfg.label         = nwaycomp.label;
      cfg.labelcolor    = 'k';
      cfg.labelsize     = 10;
      roe_brainplot_chancircle(cfg);
      title(plottitle)
      
      
      % plot B
      subplot(3,4,4)
      plot(nwaycomp.freq,B(:,icomp));
      axis tight
      ylim = get(gca,'ylim');
      set(gca,'ylim',[0 ylim(2)*1.05])
      xlabel('frequency (Hz)')
      ylabel('loading')
      title('B')
      
      % plot C
      subplot(3,4,8)
      plot(1:ntrial,C(:,icomp))
      ylim = get(gca,'ylim');
      axis tight
      set(gca,'ylim',[0 ylim(2)*1.05])
      xlabel('trial')
      ylabel('loading')
      title('C')
  
      % C t-test: rememberd vs forgotten
      subplot(3,4,12)
      ind1 = nwaycomp.trialinfo(:,3)==1; % remembered
      ind2 = nwaycomp.trialinfo(:,3)==0; % forgotten
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
      legend({['rem-forgot ' num2str(sum(ind1)) '/' num2str(sum(ind2))]})
      [h p,ci,stats12] = ttest2(C(ind1,icomp),C(ind2,icomp));
      title(['rem/forgot t = ' num2str(stats12.tstat,'%1.2f')])
    
      % set to rotation
      rotate3d on
    end
    
    
    
  end % exist 
end % isubj








