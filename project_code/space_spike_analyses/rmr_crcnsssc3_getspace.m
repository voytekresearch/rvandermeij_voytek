function rmr_crcnsssc3_getspace







% get info
info = rmr_crcnsssc3_info;

% set name suffix
fnnamesuffix = 'four_conv_timeavg_timwin20ms_freq50-50-1000';

% nway settings
nwayalg      = 'spacetime'; % 'spacefsp'
nwaynmethod  = 'splitrel'; %'ncomp15sr'; % 'splitrel'
nwaynrand    = 50; % 'rnd40_conv1e-6'
nwayconvcrit = 1e-8;
normmethod   = '8throotpower'; % 'coh'   'avgoverepoch' 'propodtod' 'kthrootpower'
nwaysplit    = '4ths'; %  oddevenspikes  oddeventrials
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];

% read data and get fourier
distcompsys = 'torque';
dosplitrel = strcmp(nwaynmethod,'splitrel') || strcmp(nwaynmethod(end-1:end),'sr');
if dosplitrel
  nwayadd = [nwayadd '_' 'split' nwaysplit];
end
for    isess = 11%:-1:1     ;% :numel(info.session); % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  disp(['working on ' currsess ', ' fnnamesuffix ', ' nwayadd])
  
  % set nway filename and continue if it doesn't exist yet
  nwayfn = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '.mat'];
  if ~exist(nwayfn,'file')
    
    % set fns and load data
    fourfn = [info.savepath currsess '_' 'fourier' '_' fnnamesuffix '.mat'];
    datafn = [info.savepath currsess '_' 'dataetc' '_' fnnamesuffix '.mat'];
    load(datafn)
    
    % create normed fourier
    fourndfn = [info.scratchpath currsess '_' 'fourier' '_' fnnamesuffix '_' normmethod '.mat'];
    load(fourfn)
    %     if strncmp(fliplr(normmethod),fliplr('throotpower'),11)
    %       [fourier,scaling] = roe_fouriernormalize(fourier,normmethod);
    %     else
    %       fourier = roe_fouriernormalize(fourier,normmethod);
    %     end
    fourier = roe_fouriernormalize(fourier,normmethod);
    save(fourndfn,'fourier','-v7.3')
    clear fourier
    
    
    % create split and save
    if dosplitrel
      switch nwaysplit
        
        case 'oddeventrials'
          % load from disk
          load(fourndfn)
          
          % set split indices (split over laps)
          lapind = data.trialinfo(:,1);
          nlap   = numel(unique(lapind)); % assume no laps are missing
          if size(fourier,3)~=nlap && (size(fourier,3)/nlap)>6
            % input is segment specific
            lapsprt1 = sort([1:4:(nlap-rem(nlap,4)) 2:4:(nlap-rem(nlap,4))]); % ensure equal ammount of laps in both sets and ensure they're approx. equal in L/R
            lapsprt2 = sort([3:4:(nlap-rem(nlap,4)) 4:4:(nlap-rem(nlap,4))]);
            if (rem(nlap,4) == 2) || (rem(nlap,4) == 3)
              lapsprt1 = [lapsprt1 nlap-rem(nlap,4)+1];
              lapsprt2 = [lapsprt2 nlap-rem(nlap,4)+2];
            end
            % get trial indices belonging to each lap
            trlindprt1 = ismember(lapind,lapsprt1);
            trlindprt2 = ismember(lapind,lapsprt2);
          elseif size(fourier,3)~=nlap && (size(fourier,3)/nlap)<=6
            % low number of segments, split over full laps
            shnlap = floor(nlap/2)*2;
            trlindprt1 = ismember(lapind,1:2:shnlap);
            trlindprt2 = ismember(lapind,2:2:shnlap);
          else
            % input is averaged over segments
            trlindprt1 = 1:2:nlap;
            trlindprt2 = 2:2:nlap;
            lowestn = min(numel(trlindprt1),numel(trlindprt2));
            trlindprt1 = trlindprt1(1:lowestn);
            trlindprt2 = trlindprt2(1:lowestn);
          end
          % save splits
          fourierfull = fourier;
          % set fns
          fourprt1fn = [info.scratchpath currsess '_' 'fourier' '_' fnnamesuffix '_' normmethod '_' nwaysplit '_' 'prt1' '.mat'];
          fourprt2fn = [info.scratchpath currsess '_' 'fourier' '_' fnnamesuffix '_' normmethod '_' nwaysplit '_' 'prt2' '.mat'];
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
          
        case 'oddevenspikes'
          % set fns
          fourprt1fn = [info.scratchpath currsess '_' 'fourier' '_' fnnamesuffix '_' normmethod '_' nwaysplit '_' 'prt1' '.mat'];
          fourprt2fn = [info.scratchpath currsess '_' 'fourier' '_' fnnamesuffix '_' normmethod '_' nwaysplit '_' 'prt2' '.mat'];
          % load base, norm, and save
          load([fourfn(1:end-4) '_' 'prunemode' '_' 'odd' '.mat']);
          %     if strncmp(fliplr(normmethod),fliplr('throotpower'),11)
          %       [fourier,scaling] = roe_fouriernormalize(fourier,normmethod);
          %     else
          %       fourier = roe_fouriernormalize(fourier,normmethod);
          %     end
          fourier = roe_fouriernormalize(fourier,normmethod);
          save(fourprt1fn,'fourier','-v7.3')
          clear fourier
          load([fourfn(1:end-4) '_' 'prunemode' '_' 'even' '.mat']);
          %     if strncmp(fliplr(normmethod),fliplr('throotpower'),11)
          %       [fourier,scaling] = roe_fouriernormalize(fourier,normmethod);
          %     else
          %       fourier = roe_fouriernormalize(fourier,normmethod);
          %     end
          fourier = roe_fouriernormalize(fourier,normmethod);
          save(fourprt2fn,'fourier','-v7.3')
          clear fourier
          % save fns for fourierdata
          fourierpart = {fourprt1fn,fourprt2fn};
          
        case '4ths'
          fourierpart = [];
          % load evert base, norm, and save
          for ipart = 1:4
            fourierpart{ipart} = [info.scratchpath currsess '_' 'fourier' '_' fnnamesuffix '_' normmethod '_' nwaysplit '_' 'prt' num2str(ipart) '.mat'];
            load([fourfn(1:end-4) '_' 'prunemode' '_' num2str(ipart) 'of4' '.mat']);
            % set scaling as normmethod if needed
            %     if strncmp(fliplr(normmethod),fliplr('throotpower'),11)
            %       [fourier,scaling] = roe_fouriernormalize(fourier,normmethod);
            %     else
            %       fourier = roe_fouriernormalize(fourier,normmethod);
            %     end
            fourier = roe_fouriernormalize(fourier,normmethod);
            save(fourierpart{ipart},'fourier','-v7.3')
            clear fourier
          end
          
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
    cfg.checkpointpath     = [info.scratchpath 'CHECKPOINT' '_' 'hc5' '_' currsess '_' nwaynmethod '_' normmethod '_'];
    % qsub options
    switch distcompsys
      case 'torque'
        cfg.distcomp.system          = 'torque';
        cfg.distcomp.timreq          = 60*60*4; %
        cfg.distcomp.matlabcmd       = '/opt/matlab/2015a/bin/matlab';
        cfg.distcomp.torquequeue     = 'hotel';
        cfg.distcomp.torquestack     = 5;
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
      cfg.ncompeststart       = 5;
      cfg.ncompeststep        = 1;
      cfg.ncompestend         = 8;
      cfg.ncompestsrcritval   = [.7 .7 0 .7 0];
      cfg.ncompestsrcritjudge = 'meanoversplitscons';
    else
      error('...')
    end
    % go for it
    nwaycomp = nd_nwaydecomposition(cfg,fourierdata);
    
    
    
    
    
    %%%%%%%%%%%%%%
    % get C's again
    if ~isempty(strfind(normmethod,'avgoverepoch'))
      
      % first get four
      normmethodC = normmethod(14:end);
      normmethodC = [normmethodC '_' 'powertrialnorm']; % add trial normalization
      load(fourfn)
      fourier = roe_fouriernormalize(fourier,normmethodC);
      
      % extract comp
      comp  = nwaycomp.randomstat.startvalglobmin;
      ncomp = numel(nwaycomp.comp);
      
      % get C
      startval = comp;
      startval{1} = bsxfun(@times,startval{1},startval{3});
      startval{3} = [];
      comp    = cell(1,nwaynrand);
      expvar  = zeros(1,nwaynrand);
      scaling = cell(1,nwaynrand);
      for irand = 1:nwaynrand
        [comp{irand},dum,dum,expvar(irand),scaling{irand}] = nwaydecomp_spacetime(fourier,ncomp,freqoi,'Dmode','identity','startval',startval,'holdparam',[0 1 0 0 1],'convcrit',nwayconvcrit,'niter',cfg.numiter,'dispprefix',['2ndpassC irand = ' num2str(irand) '/' num2str(nwaynrand) ': ']);
      end
      
      % sort and extract
      [dum, sortind] = sort(expvar,'descend');
      comp    = comp{sortind(1)};
      expvar  = expvar(sortind(1));
      scaling = scaling{sortind(1)};
      
      % put in nwaycomp and safe
      outputcomp = [];
      for icomp = 1:ncomp
        % A,B,C,S
        for iparam = 1:4
          outputcomp{icomp}{iparam} = comp{iparam}(:,icomp);
        end
        % D
        outputcomp{icomp}{5} = comp{5}(:,icomp);
      end
      comp = outputcomp;
      nwaycomp.comp    = comp;
      nwaycomp.scaling = scaling;
      % keep randomstat and splitrelstat
      
      % save
      ind = strfind(nwayfn,nwayalg);
      nwayfn = [nwayfn(1:ind-1) '2ndpassC' '_' nwayfn(ind:end)];
      save(nwayfn,'nwaycomp','-v7.3')
    end
    %%%%%%%%%%%%%%
    
    
    
    
    
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
end %isess






% get out of here
exit





















function playground




% get info
info = rmr_crcnsssc3_info;
info.datapath = '/Users/roemervandermeij/Work/Data/tmpdata/CRCNS_SSC3/';
info.savepath = '/Users/roemervandermeij/Work/Data/tmpdata/crcnsssc3/';

% set name suffix
fnnamesuffix = 'four_conv_timeavg_timwin20ms_freq50-50-1000';

% set cgname suffic
cgnamesuffix = 'contcorrgram_-20to20ms_keeptrials';

% nway settings
nwayalg      = 'spacetime'; % 'spacefsp'
nwaynmethod  = 'ncomp10sr'; %'ncomp15sh'; % 'splitrel' 'ncomp10sh'
nwaynrand    = 50; % 'rnd40_conv1e-6'
nwayconvcrit = 1e-8;
normmethod   = '8throotpower'; % % 'coh'   'avgoverepoch' 'propodtod' 'kthrootpower'
nwaysplit    = '4ths'; %  oddevenspikes  oddeventrials  4ths
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];

% read data and get fourier
if (strcmp(nwaynmethod,'splitrel') || strcmp(nwaynmethod(end-1:end),'sr'))
  nwayadd = [nwayadd '_' 'split' nwaysplit];
end
for    isess = 11%:-1:1     ;% :numel(info.session); % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  
  
  % set nway filename and continue if it doesn't exist yet
  nwayfn = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '.mat'];
  cgfn   = [info.savepath currsess '_' cgnamesuffix '.mat'];
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
    
    % compute firing rate
    % read in data
    data = rmr_crcnsssc3_readspike([info.datapath currsess],50);
    
    % select units with firing rate at least 1Hz in half of trials
    trllength = cat(1,data.time{:});
    trllength = trllength(:,2) - trllength(:,1);
    trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
    trlspikes = full(cat(2,trlspikes{:}));
    trlspikerate = trlspikes ./ trllength.';
    selind = sum(trlspikerate>1,2)>round(numel(data.trial)/2);
    for itrial = 1:numel(data.trial)
      data.trial{itrial} = data.trial{itrial}(selind,:);
    end
    data.label = data.label(selind);
  
    
    % plot
    cfg = [];
    cfg.condind   = {data.trialinfo<=round(numel(data.trial)/2),data.trialinfo>round(numel(data.trial)/2)};
    cfg.condlabel = {'1sthalf','2ndhalf'};
    cfg.cgfn      = cgfn;
    cfg.cgramdisp = 'pasquale2008'; % pasquale2008 sum
    rmr_plotspacespike(cfg,nwaycomp,data);

 
        
    
  end % exist
 
end % isess











































