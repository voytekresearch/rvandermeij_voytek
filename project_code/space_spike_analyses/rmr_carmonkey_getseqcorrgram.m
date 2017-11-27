function rmr_carmonkey_getseqcorrgram



% get info
info = rmr_carmonkey_info;
info.datapath = '/Users/roemer/Work/Data/tmpdata/carmena_monkey/';


% set name suffix
%fnnamesuffix = 'four_direct_timeavg_1hzfiring_timwin10ms_freq50-50-1000_gotilltar';
%fnnamesuffix = 'four_direct_timeavg_1hzfiring_timwin10ms_freq50-50-1000_centertillgo';
%fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin20ms_freq50-50-1000_gotilltar';
fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin20ms_freq50-50-1000_centertillgo';
fnnamesuffix = 'four_conv_timeavg_1hzfiring_timwin100ms_freq10-10-200_gotilltar';

% set nway settings
nwayalg     = 'spacetime'; % 'spacefsp'
nwaynmethod = 'splitrel'; %'ncomp15sh'; % 'splitrel' 'ncomp10sh'
nwaynrand    = 50;
nwayconvcrit = 1e-8;
normmethod   = '8throotpower'; % 'coh'   'avgoverepoch' '16throotpower'
nwaysplit    = 'oddevenspikes'; %  oddevenspikes  oddeventrials
nwayadd = [normmethod '_' nwayalg '_' nwaynmethod '_' 'rnd' num2str(nwaynrand) '_' 'conv' num2str(nwayconvcrit)];


% extract trialdeftype
ind = strfind(fnnamesuffix,'_');
trialdeftype = fnnamesuffix(ind(end)+1:end);

% read data and get fourier, 1 session
for      isubj = 1     :numel(info.subj)
  for    isess = 1     ; % 1 session for now
    didsplitrel = strcmp(nwaynmethod,'splitrel') || strcmp(nwaynmethod(end-1:end),'sr');
    if didsplitrel
      nwayadd = [nwayadd '_' 'split' nwaysplit];
    end
    
    % set currs
    currsubj = info.subj{isubj};
    currsess = info.session.(currsubj){isess};
    
    % get nway
    nwayfn = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '.mat'];
    load(nwayfn)
    
    
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
    
    
    % get sequence corrgram per component
    ncomp = numel(nwaycomp.comp);
    corrgram = cell(1,ncomp);
    for icomp = 1:ncomp
      % get units to investigate
      A = nwaycomp.comp{icomp}{1};
      [dum sortind] = sort(A,'descend');
      
      % get corrgram
      cfg = [];
      cfg.binedges   = (-20:1:20) ./ 1000; % in s
      cfg.keeptrials = 'yes';
      cfg.selectunit = data.label(sortind(1:4));
      corrgram{icomp} = rmr_seqcorrgramspiketrain(cfg,data);
    end
    
    
    % save
    data = rmfield(data,'trial');
    timebins = cfg.binedges(2:end)-mean(diff(cfg.binedges))/2;
    binedges = cfg.binedges;
    timbinname = [num2str(binedges(1).*1000) 'msto' num2str(binedges(end).*1000) 'ms' num2str(mean(diff(binedges)).*1000) 'msbins'];
    if strcmp(cfg.keeptrials,'yes')
      fn = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '_' 'seqcorrgram' '_' timbinname '_' 'keeptrials' '.mat'];
    else
      fn = [info.savepath currsess '_' fnnamesuffix '_' nwayadd '_' 'seqcorrgram' '_' timbinname  '.mat'];
    end
    save(fn,'corrgram','nwaycomp','timebins','binedges','data','-v7.3');
  end
end
















function playground

