function rmr_carmonkey_getspace


% lame
%cd('/oasis/tscc/scratch/roevdmei/');

% get info
info = rmr_carmonkey_info;
%info.datapath = '/Users/roemer/Work/Data/tmpdata/carmena_monkey/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/carmena_monkey/';


% set name suffix
%fnnamesuffix = 'four_direct_timeavg_1hzfiring_timwin10ms_freq50-50-1000_gotilltar';
%fnnamesuffix = 'four_direct_timeavg_1hzfiring_timwin10ms_freq50-50-1000_centertillgo';
%fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin20ms_freq50-50-1000_gotilltar';
%fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin20ms_freq50-50-1000_centertillgo';
%fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin100ms_freq10-10-200_gotilltar';
fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin20ms_freq50-50-1000_gotilltar';

% nway settings
nwayalg     = 'spacetime'; % 'spacefsp'
nwaynmethod = 'ncomp18'; %'ncomp15sh'; % 'splitrel' 'ncomp10sr'
nwaynrand    = 50;
nwayconvcrit = 1e-8;
normmethod   = 'avgoverepoch_8throotpower'; % 'coh'   'avgoverepoch' '8throotpower' 'kthrootpowertrialnorm' 
nwaysplit    = 'oddevenspikes'; %  oddevenspikes  oddeventrials
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];

% read data and get fourier, 1 session
distcompsys = 'torque';
dosplitrel = strcmp(nwaynmethod,'splitrel') || strcmp(nwaynmethod(end-1:end),'sr');
if dosplitrel
  nwayadd = [nwayadd '_' 'split' nwaysplit];
end
for      isubj = 1     :numel(info.subj)
  for    isess = 2     :numel(info.session.(info.subj{isubj}))
    
    % set currs
    currsubj = info.subj{isubj};
    currsess = info.session.(currsubj){isess};
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
      
      %       % process fourierdata wrt to distcomp
      %       switch distcompsys
      %         case 'torque'
      %           % nothing necessary
      %         case 'matlabpct'
      %           load(fourierdata.fourier)
      %           fourierdata.fourier = fourier;
      %           clear fourier
      %           if dosplitrel
      %             for isplit = 1:numel(fourierdata.fourierpart)
      %               load(fourierdata.fourierpart{isplit});
      %               fourierdata.fourierpart{isplit} = fourier;
      %               clear fourier
      %             end
      %           end
      %         otherwise
      %           error('woopsie')
      %       end
      
      % nway settings
      cfg = [];
      cfg.model              = nwayalg;
      cfg.datparam           = 'fourier';
      cfg.Dmode              = 'identity';
      cfg.numiter            = 3000;
      cfg.convcrit           = nwayconvcrit;
      cfg.randstart          = nwaynrand;
      cfg.outputfile         = nwayfn;
      cfg.checkpointpath     = [info.scratchpath 'CHECKPOINT' '_' 'cm' '_' currsess '_' nwaynmethod '_' normmethod '_'];
      % qsub options
      switch distcompsys
        case 'torque'
          cfg.distcomp.system          = 'torque';
          cfg.distcomp.timreq          = 60*60*2; %
          cfg.distcomp.matlabcmd       = '/opt/matlab/2015a/bin/matlab';
          cfg.distcomp.torquequeue     = 'hotel';
          cfg.distcomp.torquestack     = 10;
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
        cfg.ncompeststart       = 18;
        cfg.ncompeststep        = 6;
        cfg.ncompestend         = 20;
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
end % isubj





% get out of here
exit












function playground




% get info
info = rmr_carmonkey_info;
info.datapath = '/Users/roemervandermeij/Work/Data/tmpdata/carmena_monkey/';
info.savepath = '/Users/roemervandermeij/Work/Data/tmpdata/carmonkey/';

% set name suffix
%fnnamesuffix = 'four_direct_timeavg_1hzfiring_timwin10ms_freq50-50-1000_gotilltar';
%fnnamesuffix = 'four_direct_timeavg_1hzfiring_timwin10ms_freq50-50-1000_centertillgo';
%fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin20ms_freq50-50-1000_centertillgo';
%fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin20ms_freq50-50-1000_centertillgo';
fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin100ms_freq10-10-200_gotilltar';
fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin20ms_freq50-50-1000_gotilltar';

% set cgname suffic
cgnamesuffix = ['1hzfiring_gotilltar' '_' 'contcorrgram' '_' '-20to20ms_keeptrials'];


% extract trialdeftype
ind = strfind(fnnamesuffix,'_');
trialdeftype = fnnamesuffix(ind(end)+1:end);

% nway settings
nwayalg     = 'spacetime'; % 'spacefsp'
nwaynmethod = 'splitrel'; %'ncomp15sr'; % 'splitrel' 'ncomp10sr'
nwaynrand    = 50;
nwayconvcrit = 1e-8;
normmethod   = '8throotpower'; % 'coh'   'avgoverepoch' '16throotpower'  powertrialnorm_8throotpower
nwaysplit    = 'oddevenspikes'; %  oddevenspikes  oddeventrials 4ths
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];

% read data and get fourier, 1 session
if (strcmp(nwaynmethod,'splitrel') || strcmp(nwaynmethod(end-1:end),'sr'))
  nwayadd = [nwayadd '_' 'split' nwaysplit];
end
for      isubj = 1     :numel(info.subj)
  for    isess = 1     :numel(info.session.(info.subj{isubj}))
 
    % set currs
    currsubj = info.subj{isubj};
    currsess = info.session.(currsubj){isess};
     
    
    % set nway filename and continue if it doesn't exist yet
    nwayfn = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '.mat'];
    cgfn   = [info.savepath currsess '_' cgnamesuffix '.mat'];
    if exist(nwayfn,'file')
      load(nwayfn)
      
      
      % plot decomp stats
      cfg = [];
      cfg.savefig    = 'no';
      cfg.figprefix  = 'paco';
      cfg.figvisible = 'on';
      roe_nwaystatplot(cfg,nwaycomp)
       
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
       
      % fetch data
      fsample = 40000; % oversample intentionally
      data = rmr_carmonkey_readspike([info.datapath currsess],fsample,trialdeftype);
      % select units with firing rate above 1 hz
      trllength = cat(1,data.time{:});
      trllength = trllength(:,2) - trllength(:,1);
      trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
      trlspikes = full(cat(2,trlspikes{:}));
      spikerate = sum(trlspikes,2) ./ sum(trllength);
      selind = spikerate>1;
      for itrial = 1:numel(data.trial)
        data.trial{itrial} = data.trial{itrial}(selind,:);
      end
      data.label = data.label(selind);
       
      % plot
      cfg = [];
      cfg.condind   = {nwaycomp.trialinfo(:,1)~=12,nwaycomp.trialinfo(:,1)==12}; % 7 == entering target, 12 == timeout
      cfg.condlabel = {'cursor reached target (hit)','failed to reach target (miss)'};
      cfg.cgfn      = cgfn;
      cfg.cgramdisp = 'sum'; % pasquale2008 sum
      rmr_plotspacespike(cfg,nwaycomp,data);
      
      
       
    end % exist
     
     
  end % isess
end






































% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % set cgname suffic
% cgnamesuffix = 'corrgram_1hzfiring_gotilltar_-20msto20ms1msbins_keeptrials';
% 
% % read data and get fourier, 1 session
% for      isubj = 1     :numel(info.subj)
%   for    isess = 1     ; % 1 session for now
%     
%     % set nway filename and continue if it doesn't exist yet
%     nwayfn = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '.mat'];
%     cgfn   = [info.savepath currsess '_' cgnamesuffix '.mat'];
%     if exist(nwayfn,'file')
%       load(nwayfn)
%       load(cgfn)
%       
%       % fetch data
%       fsample = 40000; % oversample intentionally
%       data = rmr_carmonkey_readspike([info.datapath currsess],fsample,trialdeftype);
%       % select units with firing rate above 1 hz
%       trllength = cat(1,data.time{:});
%       trllength = trllength(:,2) - trllength(:,1);
%       trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
%       trlspikes = full(cat(2,trlspikes{:}));
%       spikerate = sum(trlspikes,2) ./ sum(trllength);
%       selind = spikerate>1;
%       for itrial = 1:numel(data.trial)
%         data.trial{itrial} = data.trial{itrial}(selind,:);
%       end
%       data.label = data.label(selind);
%       trlspikes = cellfun(@sum,data.trial,repmat({2},[1 numel(data.trial)]),'uniformoutput',0);
%       trlspikes = full(cat(2,trlspikes{:}));
%       trialtime = cellfun(@max,data.time);
%       spikerate = bsxfun(@rdivide,trlspikes,trialtime);
%       
%       % set basics
%       comp   = nwaycomp.comp;
%       ncomp  = numel(comp);
%       nchan  = size(comp{1}{1},1);
%       nfreq  = size(comp{1}{2},1);
%       ntrial = size(comp{1}{3},1);
%       
%       
%       % gather A,B,C,S
%       A = zeros(nchan,ncomp);
%       B = zeros(nfreq,ncomp);
%       C = zeros(ntrial,ncomp);
%       S = zeros(nchan,ncomp);
%       for icomp = 1:ncomp
%         A(:,icomp) = comp{icomp}{1};
%         B(:,icomp) = comp{icomp}{2};
%         C(:,icomp) = comp{icomp}{3} .^ 2; %%% SQUARE
%         S(:,icomp) = comp{icomp}{4};
%       end
%       % sort by non-unitary spatial maps
%       tmpA = A ./ repmat(max(A,[],1),[size(A,1) 1]);
%       tmpA = sort(tmpA,1);
%       [tmpA comporder] = sort(tmpA(end-1,:),'descend');
%       %comporder(tmpA<0.1) = [];
%       %comporder = 1:ncomp;
%       %comporder = comporder(1:2);
%       
%             
%       % obtain STTCs
%       cfg = [];
%       cfg.binedges = [timebins-(mean(diff(timebins))/2) timebins(end)+(mean(diff(timebins))/2)];
%       cfg.trials   = {find(nwaycomp.trialinfo(:,3)~=12),find(nwaycomp.trialinfo(:,3)==12)};
%       sttc = rmr_sttcspiketrain(cfg,data,corrgram);
%    
%       % per comp
%       for icomp = comporder
%         
%         figure('numbertitle','off','name',['comp' num2str(icomp)]);
%         
%         
%         %         % plot A
%         %         subplot(4,4,1)
%         %         plot(1:nchan,A(:,icomp));
%         %         ylim = get(gca,'ylim');
%         %         ylim = round([0 ylim(2)*1.05]*10)./10;
%         %         set(gca,'ylim',ylim,'xlim',[1 nchan],'ytick',ylim)
%         %         xlabel('units')
%         %         ylabel('loading (au)')
%         %         title('spatial map of network')
%         
%         % plot A
%         subplot(4,4,1)
%         [dum, sortind] = sort(mean(spikerate,2),'descend');
%         %sortind = 1:numel(sortind);
%         [ax h1 h1] = plotyy(1:nchan,A(sortind,icomp),1:nchan,mean(spikerate(sortind,:),2));
%         ylim = get(ax(1),'ylim');
%         set(ax(1),'ylim',[0 ylim(2)*1.05],'xlim',[1 nchan])
%         set(ax(2),'ylim',[0 20],'ytick',[0:2.5:20],'xlim',[1 nchan])
%         xlabel('sorted units')
%         ylabel('loading')
%         %legend([repmat('comp',[ncomp 1]) deblank(num2str((1:ncomp)'))])
%         title('A')
%         
%         % plot C t-test: enter target or not
%         subplot(4,4,2)
%         condind1 = nwaycomp.trialinfo(:,1)~=12; % 7 == entering target, 12 == timeout
%         condind2 = ~condind1;
%         m1 = mean(C(condind1,icomp));
%         m2 = mean(C(condind2,icomp));
%         sem1 = std(C(condind1,icomp)) ./ sqrt(sum(condind1));
%         sem2 = std(C(condind2,icomp)) ./ sqrt(sum(condind2));
%         errorbar([m1 m2],[sem1 sem2])
%         ylim = get(gca,'ylim');
%         ylim = round([0 ylim(2)*1.1]*100)./100;
%         set(gca,'ylim',ylim);
%         set(gca,'xlim',[0.5 2.5])
%         set(gca,'xtick',[1 2])
%         set(gca,'ytick',ylim)
%         set(gca,'xticklabel',{'cursor reached target (hit)','failed to reach target (miss)'})
%         ylabel('loading (au)')
%         [h p,ci,stats] = ttest2(C(condind1,icomp),C(condind2,icomp));
%         title(['network strength hits/miss trials, mean+sem, t = ' num2str(stats.tstat,'%1.2f')])
%         
%         % plot S
%         subplot(4,4,3)
%         hold on
%         [dum sortind] = sort(A(:,icomp),'descend');
%         [ax h1 h2] = plotyy(1:5,(S(sortind(1:5),icomp)-S(sortind(1),icomp))*1000,1.5:4.5,diff((S(sortind(1:5),icomp)-S(sortind(1),icomp))*1000));
%         set(ax(1),'ylim',[-.010 .010]*1000,'xlim',[1 5],'xtick',1:5)
%         set(ax(1),'ytick',(-.01 :0.005: .01)*1000)
%         set(ax(2),'ylim',[-.005 .005]*1000,'xlim',[1 5],'xtick',1:5)
%         set(ax(2),'ytick',(-.005 :0.0025: .005)*1000)
%         set(h2,'linestyle','--','marker','o','markerfacecolor', [0.85 0.325 0.098])
%         xlabel('units sorted by spatial map strength')
%         ylabel(ax(1),'time delay (ms)')
%         ylabel(ax(2),'diff of time delay (ms)')
%         title('time delays of units')
%                 
%         % plot pairs of corrgrams
%         [dum sortind] = sort(A(:,icomp),'descend');
%         chanind = sortind(1:4);
%         count = 4;
%         for iseed = 1:numel(chanind)
%           for itarg = iseed:numel(chanind)
%             if iseed == 3 && itarg == 3
%               count = count+2;
%             else
%               count = count+1;
%             end
%             subplot(4,4,count);
%             setcol = {[0 .75 1],[0 0 .75]};
%             hold on
%             for iset = 1:2
%               plot(timebins*1000,squeeze(sttc(iset,chanind(iseed),chanind(itarg),:)),'linewidth',1,'color',setcol{iset});
%             end
%             set(gca,'xtick',[-20 -10 0 10 20])
%             xlabel('time from seed unit spikes (ms)')
%             ylabel('STTC')
%             title(['unit' num2str(iseed) '-' 'unit' num2str(itarg)])
%             set(gca,'ylim',[-.1 .1],'ytick',[-.1 0 .1]);
%             legend('cond1','cond2')
%           end
%         end
%         
%         
%                 
%       end % comporder
%       
%       
%       
%       
%     end % exist
%   end % isess
% end

































