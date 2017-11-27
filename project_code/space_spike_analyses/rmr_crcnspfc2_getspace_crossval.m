function rmr_crcnspfc2_getspace_crossval




% get info
info = rmr_crcnspfc2_info;
info.datapath = '/Users/roemer/Work/Data/tmpdata/CRCNS_PFC2/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/crcnspfc2/';

% set name suffix
fnnamesuffix = 'four_saw_timeavg_1segsel_timwin25ms_freq20-20-500';


% nway settings
nwayalg      = 'spacetime'; % 'spacefsp'
nwaynmethod  = 'ncomp10sr'; %'ncomp15sh'; % 'splithalf' 'ncomp10sh'
nwaynrand    = 50; % 'rnd40_conv1e-6'
nwayconvcrit = 1e-8;
normmethod   = '16throotpower'; % 'coh'   'avgoverepoch' 'propodtod' 'kthrootpower'
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];
nwaysplit    = '4ths'; %  oddevenspikes  oddeventrials
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];

% prune settings
nprune    = 4;
pruneperc = 25;

% read data and get fourier
for    isess = 2;%:-1:1     ;% :numel(info.session); % 1 session for now
  dosplitrel = strcmp(nwaynmethod,'splitrel') || strcmp(nwaynmethod(end-1:end),'sr');
  if dosplitrel
    nwayadd = [nwayadd '_' 'split' nwaysplit];
  end
    % set currs
  currsess = info.session{isess};
  
  % set nway filename and continue if it doesn't exist yet
  nwayfn = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '.mat'];
  if exist(nwayfn,'file')
    load(nwayfn)
    % do the same for fncv
    %fncv = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '_' 'crossvalnprune' num2str(nprune) '_' num2str(pruneperc) 'perc' '.mat'];
    fncv = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '_' 'crossvalnprune' num2str(nprune) '_' 'nof4' '.mat'];
    if ~exist(fncv,'file')
      
      % set basics
      ncomp  = numel(nwaycomp.comp);
      freqoi = nwaycomp.freq;
      
      % extract component loadings in low level format
      maincomp = nwaycomp.randomstat.startvalglobmin;
      
      % load original fourier to get scaling coefficients
      fourfn = [info.savepath currsess '_' 'fourier' '_' fnnamesuffix '.mat'];
      load(fourfn)
      [fourier,fourscaling] = roe_fouriernormalize(fourier,normmethod);
      
      % compute crossvalidation loadings for ABCS of every pruneset
      compcv = cell(1,nprune);
      for iprune = 1:nprune
        % load pruned fourier
        %prfn = [info.savepath currsess '_' 'fourier' '_' fnnamesuffix '_' 'pruneset' num2str(iprune) '_' num2str(pruneperc) '.mat'];
        prfn = [info.savepath currsess '_' 'fourier' '_' fnnamesuffix '_' 'prunemode' '_' num2str(iprune) 'of4' '.mat'];
        load(prfn);
        fourier = roe_fouriernormalize(fourier,fourscaling);
        
        % obtain comp
        [comp,startval,ssqres,expvar,scaling] = nwaydecomp_spacetime(fourier,ncomp,freqoi,'startval',maincomp,'Dmode','identity','convcrit',1e-8);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % match comp with maincomp
        nparam = 5;
        model  = 'spacetime';
        freq   = freqoi;
        Dmode  = 'identity';
        % set current estcomp
        currcomp = cell(1,2);
        currcomp{1} = maincomp;
        currcomp{2} = comp;
        % compute component congruence for all possible pairs between selected rands
        compcongr = zeros(ncomp,ncomp,length(currcomp{1}));
        for icompr1 = 1:ncomp
          for icompr2 = 1:ncomp
            for iparam = 1:nparam
              % perform model specific stuff
              switch model
                case {'spacetime','spacefsp'}
                  switch iparam
                    case {1,2}
                      paramc1 = currcomp{1}{iparam}(:,icompr1);
                      paramc2 = currcomp{2}{iparam}(:,icompr2);
                      % normalize
                      paramc1 = paramc1 ./ sqrt(sum(abs(paramc1).^2));
                      paramc2 = paramc2 ./ sqrt(sum(abs(paramc2).^2));
                      % put in compsh
                      compcongr(icompr1,icompr2,iparam) = abs(paramc1' * paramc2);
                    case 3
                      paramc1 = currcomp{1}{iparam}(:,icompr1);
                      paramc2 = currcomp{2}{iparam}(:,icompr2);
                      % normalize
                      paramc1 = paramc1 ./ sqrt(sum(abs(paramc1).^2));
                      paramc2 = paramc2 ./ sqrt(sum(abs(paramc2).^2));
                      % put in compsh
                      if numel(paramc1) == numel(paramc2)
                        compcongr(icompr1,icompr2,iparam) = abs(paramc1' * paramc2);
                      else
                        compcongr(icompr1,icompr2,iparam) = 0; % congruence can't be computed, set to maximally incongruent (0)
                      end
                    case 4
                      switch model
                        case 'spacetime'
                          % create frequency specific phases weighted by spatial maps and frequency profiles
                          A1 = currcomp{1}{1}(:,icompr1);
                          A2 = currcomp{2}{1}(:,icompr2);
                          B1 = currcomp{1}{2}(:,icompr1);
                          B2 = currcomp{2}{2}(:,icompr2);
                          S1 = currcomp{1}{4}(:,icompr1);
                          S2 = currcomp{2}{4}(:,icompr2);
                          % construct complex site by freq matrix
                          Scomp1 = exp(1i*2*pi*repmat(freq(:).',[size(A1,1) 1]).*repmat(S1,[1 size(B1,1)]));
                          Scomp2 = exp(1i*2*pi*repmat(freq(:).',[size(A2,1) 1]).*repmat(S2,[1 size(B2,1)]));
                          % scale with A
                          Scomp1 = Scomp1 .* repmat(A1,[1 size(B1,1)]);
                          Scomp2 = Scomp2 .* repmat(A2,[1 size(B2,1)]);
                          % compute congruence over freqs, than abs, then average weighted with B
                          shoverfreq = zeros(numel(B1),1);
                          for ifreq = 1:numel(B1)
                            currS1 = Scomp1(:,ifreq);
                            currS2 = Scomp2(:,ifreq);
                            currS1 = currS1 ./ sqrt(sum(abs(currS1).^2)); % not necessary now, but just in case we ever decide to not-normalize A
                            currS2 = currS2 ./ sqrt(sum(abs(currS2).^2));
                            shoverfreq(ifreq) = abs(currS1'*currS2);
                          end
                          shsumfreq = sum(shoverfreq .* (B1.*B2)) ./ sum(B1.*B2);
                          % put in compsh
                          compcongr(icompr1,icompr2,iparam) = shsumfreq;
                      end
                    case 5
                      switch Dmode
                        case 'identity'
                          % D is fixed with arbitrary order, make its congruence coefficient irrelevant
                          compcongr(icompr1,icompr2,iparam) = 1;
                      end
                  end
              end
            end
          end
        end
        % get cong coefficients by selecting most-similair unique pairings
        compcongrsel = zeros(ncomp,nparam);
        congrsum     = sum(compcongr,3);
        % match from perspective of first rand (i.e. find components of rand 2 that match those of rand 1)
        % do so by starting from the component-pair with the highest similarity, then the next most similar, etc.
        r1ind   = zeros(1,ncomp);
        r2ind   = zeros(1,ncomp);
        for icomp = 1:ncomp
          [dum, r1ind(icomp)] = max(max(congrsum,[],2));
          [dum, r2ind(icomp)] = max(congrsum(r1ind(icomp),:));
          congrsum(r1ind(icomp),:) = 0;
          congrsum(:,r2ind(icomp)) = 0;
        end
        % sanity check
        if any(diff(sort(r1ind))==0) || any(diff(sort(r2ind))==0)
          error('some components were selected multiple times')
        end
        % sort
        [r1ind, sortind] = sort(r1ind);
        r2ind = r2ind(sortind);
        
        % r2ind is the ind for comp
        comp{1} = comp{1}(:,r2ind);
        comp{2} = comp{2}(:,r2ind);
        comp{3} = comp{3}(:,r2ind);
        comp{4} = comp{4}(:,r2ind);
        comp{5} = comp{5}(:,r2ind);
        
        % store the crossval loadings
        compcv{iprune} = comp;
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % prep C and S for computing congruence coefficient
      % norm C
      maincomp{3} = maincomp{3} ./ repmat(sqrt(sum(abs(maincomp{3}).^2,1)),[size(maincomp{3},1) 1]);
      % convert S to complex domain and weight with A
      sigmacirc = 0.050;
      radconv   = ((2*pi)./sigmacirc);
      Scomp = exp(1i*maincomp{4}.*radconv);
      Sweighting = maincomp{1} ./ repmat(sum(maincomp{1},1),[size(maincomp{1},1) 1]);
      
      % mean center and renorm ABC
      for iparam = 1:3
        maincomp{iparam} = maincomp{iparam} - repmat(mean(maincomp{iparam},1),[size(maincomp{iparam},1) 1]);
        maincomp{iparam} = maincomp{iparam} ./ repmat(sqrt(sum(abs(maincomp{iparam}).^2,1)),[size(maincomp{iparam},1) 1]);
      end
      for iprune = 1:nprune
        for iparam = 1:3
          compcv{iprune}{iparam} = compcv{iprune}{iparam} - repmat(mean(compcv{iprune}{iparam},1),[size(compcv{iprune}{iparam},1) 1]);
          compcv{iprune}{iparam} = compcv{iprune}{iparam} ./ repmat(sqrt(sum(abs(compcv{iprune}{iparam}).^2,1)),[size(compcv{iprune}{iparam},1) 1]);
        end
      end
      
      % then, compute a congruence coefficient for every crossval parameter
      cvcongr  = cell(1,4);
      cvcongr{1} = NaN(ncomp,nprune);
      cvcongr{2} = NaN(ncomp,nprune);
      cvcongr{3} = NaN(ncomp,nprune);
      cvcongr{4} = NaN(ncomp,nprune);
      for iprune = 1:nprune
        % A
        cvcongr{1}(:,iprune) = diag(maincomp{1}'*compcv{iprune}{1});
        % A
        cvcongr{2}(:,iprune) = diag(maincomp{2}'*compcv{iprune}{2});
        % C
        cvcongr{3}(:,iprune) = diag(maincomp{3}'*compcv{iprune}{3});
        % A
        cvcongr{4}(:,iprune) = abs(sum((Scomp .* conj(exp(1i*compcv{iprune}{4}.*radconv))) .* Sweighting,1));
      end
      
      
      % save that shit
      save(fncv,'cvcongr','compcv','-v7.3')
    end % exist
    
  end % exist
end % isess

















function playground




% get info
info = rmr_crcnspfc2_info;
info.datapath = '/Users/roemer/Work/Data/tmpdata/CRCNS_PFC2/';
%info.savepath = '/Users/roemer/Work/Data/tmpdata/crcnspfc2/';

% set name suffix
fnnamesuffix = 'four_saw_timeavg_1segsel_timwin25ms_freq20-20-500';


% nway settings
nwayalg      = 'spacetime'; % 'spacefsp'
nwaynmethod  = 'ncomp10'; %'ncomp15sh'; % 'splithalf' 'ncomp10sh'
nwaynrand    = 40; % 'rnd40_conv1e-6'
nwayconvcrit = 1e-10;
normmethod   = '16throotpower'; % 'coh'   'avgoverepoch'
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];

% prune settings
nprune    = 20;
pruneperc = 10;

% read data and get fourier
for    isess = 2%:-1:1     ;% :numel(info.session); % 1 session for now
  
  % set currs
  currsess = info.session{isess};
  
  % set nway filename ingredients
  nwaynmethod  = {'ncomp10', 'ncomp20', 'ncomp30'};
  nwaynrand    = [40 40 40];
  
  % load
  for inway = 1:numel(nwaynmethod)
    nwayadd = [normmethod '_' nwayalg '_' nwaynmethod{inway} '_' 'rnd' num2str(nwaynrand(inway)) '_' 'conv' num2str(nwayconvcrit)];
    nwayfn = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '.mat'];
    load(nwayfn)
    fncv   = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '_' 'crossvalnprune' num2str(nprune) '_' num2str(pruneperc) 'perc' '.mat'];
    load(fncv)

    
    %
    ncomp   = numel(nwaycomp.comp);
    plotname = {'A','B','C','S'};
    figure('numbertitle','off','name',[currsess '_' 'ncomp' num2str(ncomp)])
    for iparam = 1:4
      subplot(2,2,iparam)
      imagesc(cvcongr{iparam})
      axis square
      if iparam~=4
        caxis([-1 1])
      else
        caxis([0 1])
      end
      xlabel('pruned dataset')
      ylabel('component')
      colorbar
      title(plotname{iparam})
    end
  end
end



















