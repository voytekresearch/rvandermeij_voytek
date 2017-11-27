function rmr_abcaudecog_getpsdchifreqwphase_sh_eqntrl_decomp


% get details
info = rmr_abcaudecog_info;

% set stuff
psdpath   = [info.savepath 'psd/'];

% nway settings
nwaynmethod = 'splithalf';
nwayalg     = 'parafac';


% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  % skip subj if done already
  fnfreqwphase = [psdpath currsubj '_' 'att-ign' '_freqwphase_sh_eqntrl_3cyc75msoffschi80-500.mat'];
  %fnfreqwphase = [psdpath currsubj '_' 'att-ign' '_freqwphase_sh_eqntrl_3cyc75msoffschi.mat'];
  if ~exist(fnfreqwphase,'file')
    continue
  else
    load(fnfreqwphase);
  end
 
  % set nway filename and check for precense
  nwayfn = [psdpath currsubj '_' 'freqwphase_sh_eqntrl_3cyc75msoffschi80-500' '_' nwayalg '_' 'nwaycomp' '_' nwaynmethod '.mat'];
  if exist(nwayfn,'file')
    %continue
  end
  

  
  % switch over possibilities
  cfg = [];
  cfg.numitt            = 1000;
  cfg.convcrit          = 1e-8;
  cfg.randstart         = 20;
  cfg.ncompestrandstart = 20;
  cfg.model             = nwayalg;
  cfg.complexdims       = [1 1 0 0];
  cfg.datparam          = 'wplf';

%   % select WITHIN CHAN PAC ONLY
%   nchan = numel(freqwphase.label_old);
%   newwplfsize = size(freqwphase.wphaselockspctrm);
%   newwplfsize = newwplfsize([2 3 4]);
%   freqwphase.newwphaselockspctrm = NaN(newwplfsize);
%   freqwphasesh1.newwphaselockspctrm = NaN(newwplfsize);
%   freqwphasesh2.newwphaselockspctrm = NaN(newwplfsize);
%   for ichan = 1:nchan
%     freqwphase.newwphaselockspctrm(ichan,:,:)    = squeeze(freqwphase.wphaselockspctrm(ichan,ichan,:,:));
%     freqwphasesh1.newwphaselockspctrm(ichan,:,:) = squeeze(freqwphasesh1.wphaselockspctrm(ichan,ichan,:,:));
%     freqwphasesh2.newwphaselockspctrm(ichan,:,:) = squeeze(freqwphasesh2.wphaselockspctrm(ichan,ichan,:,:));    
%   end
%   freqwphase.wphaselockspctrm    = freqwphase.newwphaselockspctrm;
%   freqwphasesh1.wphaselockspctrm = freqwphasesh1.newwphaselockspctrm;
%   freqwphasesh2.wphaselockspctrm = freqwphasesh2.newwphaselockspctrm;
%   freqwphase    = rmfield(freqwphase,'newwphaselockspctrm');
%   freqwphasesh1 = rmfield(freqwphasesh1,'newwphaselockspctrm');
%   freqwphasesh2 = rmfield(freqwphasesh2,'newwphaselockspctrm');
  
  
  % build fourierdata
  freqdata = [];
  freqdata.wplf        = freqwphase.wphaselockspctrm;
  freqdata.wplfpart{1} = freqwphasesh1.wphaselockspctrm;
  freqdata.wplfpart{2} = freqwphasesh2.wphaselockspctrm;
  freqdata.phasfreq    = freqwphase.phasfreq;
  freqdata.label       = freqwphase.label_old;
  freqdata.dimord      = freqwphase.dimord;
  freqdata.cfg         = freqwphase.cfg;
  clear freqwphase freqwphasesh1 freqwphasesh2
    
  % start nwaycomp
  cfg.ncompest           = nwaynmethod;
  cfg.ncompestshcritval  = [.6 .6 0 0];
  %cfg.ncompestcorconval  = .7;
  %cfg.ncompestvarinc     = 2.5;
  cfg.ncompeststart      = 1;
  cfg.ncompestend        = 15;
  cfg.ncompeststep       = 2;
  cfg.ncompestshdatparam = 'wplfpart';
  nwaycomp = nd_nwaydecomposition(cfg,freqdata);
  
  % save decomp data
  save(nwayfn,'nwaycomp')
  
  % clear remaining stuff
  clear fourier nwaycomp fourierdata
   
  
  
end























function playground





% get details
info = rmr_abcaudecog_info;

% set stuff
psdpath   = [info.savepath 'psd/'];

% nway settings
nwaynmethod = 'splithalf';
nwayalg     = 'parafac';

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(['working on ' currsubj])
  
  % set nway filename and check for precense
  nwayfn = [psdpath currsubj '_' 'freqwphase_sh_eqntrl_3cyc75msoffschi' '_' nwayalg '_' 'nwaycomp' '_' nwaynmethod '.mat'];
  if ~exist(nwayfn,'file')
    continue
  else
    load(nwayfn)
  end
  
  % plot decomp stats
  cfg = [];
  cfg.savefig    = 'no';
  cfg.figprefix  = currsubj;
  cfg.figvisible = 'on';
  roe_nwaystatplot(cfg,nwaycomp)
  continue
  
  % prep plotting
  comp  = nwaycomp.comp;
  freq  = 8:30;%nwaycomp.freq; % woopsie
  label = nwaycomp.label;
  ncomp = numel(comp); 
  nchan = numel(label);
  nfreq = freq;
  load([info.savepath currsubj '_fieldtrip_layout.mat']);
  [dum chanind] = intersect(lay.label,label);
     
  
  % plot ampchan
  plotlab = {'ampchan','phaschan'};
  for iplot = 1:2
    figure('numbertitle','off','name',[currsubj ' ' plotlab{iplot}])
    for icomp = 1:ncomp
      subplot(ceil(ncomp.^.5),ceil(ncomp.^.5),icomp)
      hold on
      % plot outlines
      for ioutline = 1:numel(lay.outline)
        xline = lay.outline{ioutline}(:,1);
        yline = lay.outline{ioutline}(:,2);
        line(xline,yline,'linewidth',1,'linestyle','--')
      end
      % plot circles of A
      x = lay.pos(chanind,1);
      y = lay.pos(chanind,2);
      z = zeros(size(x));
      chanmag  = abs(comp{icomp}{iplot});
      chanphas = angle(comp{icomp}{iplot});
      chansiz = 50*chanmag;
      cfg = [];
      cfg.markercolor   = 'cparam';
      cfg.viewpoint     = [0 90];
      cfg.coordinates   = [x y z];
      cfg.cparam        = chanphas;
      cfg.renderbrain   = 'no';
      cfg.colormap      = 'hsv';
      cfg.clim          = [-pi pi];
      cfg.colorbar      = 'yes';
      cfg.axissqueeze   = 'yes';
      cfg.markersize    = num2cell(chansiz);
      roe_brainplot_chancircle(cfg);
    end
  end
  
  % plot freqprof
  figure('numbertitle','off','name',[currsubj ' freqprof'])
  for icomp = 1:ncomp
    subplot(ceil(ncomp.^.5),ceil(ncomp.^.5),icomp)
    hold on
    plot(freq,comp{icomp}{3})   
    axis tight
    xlabel('Hz')
    ylabel('loading')
  end

  % plot condprof
  figure('numbertitle','off','name',[currsubj ' condprof'])
  for icomp = 1:ncomp
    subplot(ceil(ncomp.^.5),ceil(ncomp.^.5),icomp)
    if isreal(comp{icomp}{4})
      plot([1 2],comp{icomp}{4})
    else
      plot([1 2],abs(comp{icomp}{4}))
    end
    ylabel('loading')
    set(gca,'xtick',[1 2],'xticklabel',{'att','ign'})
    ylim = get(gca,'ylim');
    set(gca,'ylim',[min([0 ylim(1)]) ylim(2)])
  end

  
  
  
end





close all
isubj = isubj +1

















































