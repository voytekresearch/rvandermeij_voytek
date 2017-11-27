function rmr_predfaceval_behav



% get details
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 6    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(' ')
  disp(['patient ' currsubj])
  if isempty(info.sessionevents.(currsubj))
    disp('data not yet present')
    continue
  end
  
  
  % gather behavior
  predunpred = [];
  fearneut   = [];
  accuracy   = [];
  rt         = [];
  ntrial     = 0;
  ntrialall  = 0;
  ntrialrej  = 0;
  nsess = numel(info.sessionevents.(currsubj));
  for    isess = 1   :nsess
    
    % set
    currevents = info.sessionevents.(currsubj){isess};
    currdata   = info.sessiondata.(currsubj){isess};
    % fetch trl
    cfg = [];
    cfg.datafile  = [info.datapath currsubj '/' currdata];
    cfg.eventfile = [info.datapath currsubj '/' currevents];
    cfg.diodechan = info.sessiondiode.(currsubj){isess};
    cfg.prestim   = 0.500; % the period, in seconds, before CUE ONSET that is additionally cut out (t=0 will remain CUE ONSET)
    cfg.poststim  = 0.500; % the period, in seconds, after FACE ONSET that is additionally cut out (FACE ONSET is kept in data.trialinfo, see above)
    cfg.debugflg  = false;
    [trlall,event] = rmr_predfaceval_definetrials(cfg); % obtain the trl matrix, which contains the segmentation details
    % reject using a trick and ft_rejectartifact
    cfg = [];
    cfg.trl = trlall;
    [dum, main, ext] = fileparts(currdata);
    cfg.headerfile = [info.datapath currsubj '/' currdata];
    cfg.artfctdef.roemer.artifact = info.(currsubj).artfctdef.visual.(([ext(2:end) main]));
    cfgout = ft_rejectartifact(cfg);
    trl    = cfgout.trl;
    
    % check if its empty :(
    if isempty(trl)
      warning([currdata ' contains no trials after artifact rejection'])
      continue
    end
    
    % The trl matrix describes a trial per row, where each column describes:
    %  1) begsample   - beginning sample of trial
    %  2) endsample   - end sample of trial
    %  3) offset      - offset in samples of t=0 from start of trial
    %  4) pred/unpred - predictive (1) or unpredictive (2) cue trial
    %  5) fear/neut   - fearful (1) or neutral (2) face trial
    %  6) hit/miss    - detected valence correctly (1) or incorrectly (0)
    %  7) RT          - reaction time in seconds
    %  8) faceonset   - onset of face stimulus in seconds (from cue onset, t=0)
    %  9) respcueons  - estimated onset of response cue in seconds (from cue onset, t=0)
    % 10) faceindex   - index of shown face as defined in PTB task script
    % 11) cuediodur   - duration of cue on the screen as measured by the photo diode
    % 12) facediodur  - duration of face+isi on the screen as measured by the photo diode
    
    % remove trials with cuefaceinterval above 1 second (next step is 1.5)
    remind = (trl(:,8) - trl(:,11)) > 1.25;
    trl(remind,:) = [];
    
    predunpred = [predunpred;  trl(:,4)];
    fearneut   = [fearneut;    trl(:,5)];
    accuracy   = [accuracy;    trl(:,6)];
    rt         = [rt;          trl(:,7)];
    ntrial     = ntrial + size(trl,1);
    ntrialrej  = ntrialrej + (size(trlall,1) - size(trl,1) - sum(remind));
    ntrialall  = ntrialall + numel(remind);
  end
  rtmed = median(rt);
  rtmad = mad(rt,1);
  remind = (rt < (rtmed-6*rtmad)) | (rt > (rtmed+6*rtmad));
  rt(remind)         = [];
  accuracy(remind)   = [];
  fearneut(remind)   = [];
  predunpred(remind) = [];
  ntrial = ntrial - sum(remind);
  
  % compute/select
  ind = predunpred==1 & fearneut == 1;
  fearpredacc   = sum(accuracy(ind)) ./ sum(ind);
  fearpredrtm   = sum(rt(ind)) ./ sum(ind);
  fearpredrts   = std(rt(ind));
  fearpredn     = sum(ind);
  %
  ind = predunpred==2 & fearneut == 1;
  fearunpredacc = sum(accuracy(ind)) ./ sum(ind);
  fearunpredrtm = sum(rt(ind)) ./ sum(ind);
  fearunpredrts = std(rt(ind));
  fearunpredn   = sum(ind);
  %
  ind = predunpred==1 & fearneut == 2;
  neutpredacc   = sum(accuracy(ind)) ./ sum(ind);
  neutpredrtm   = sum(rt(ind)) ./ sum(ind);
  neutpredrts   = std(rt(ind));
  neutpredn     = sum(ind);
  %
  ind = predunpred==2 & fearneut == 2;
  neutunpredacc = sum(accuracy(ind)) ./ sum(ind);
  neutunpredrtm = sum(rt(ind)) ./ sum(ind);
  neutunpredrts = std(rt(ind));
  neutunpredn   = sum(ind);
  %
  ind = predunpred==1;
  predacc       = sum(accuracy(ind)) ./ sum(ind);
  predrtm       = sum(rt(ind)) ./ sum(ind);
  predrts       = std(rt(ind));
  %
  ind = predunpred==2;
  unpredacc     = sum(accuracy(ind)) ./ sum(ind);
  unpredrtm     = sum(rt(ind)) ./ sum(ind);
  unpredrts     = std(rt(ind));
  %
  ind = fearneut==1;
  fearacc       = sum(accuracy(ind)) ./ sum(ind);
  fearrtm       = sum(rt(ind)) ./ sum(ind);
  fearrts       = std(rt(ind));
  %
  ind = fearneut==2;
  neutacc       = sum(accuracy(ind)) ./ sum(ind);
  neutrtm       = sum(rt(ind)) ./ sum(ind);
  neutrts       = std(rt(ind));
  %
  % d prime = norminv(HR) - norminv(FAR)
  % c prime =  -(norminv(HR) + norminv(FAR))/2
  
  % spit it out
  disp([currsubj ': ' num2str(ntrial) ' usable trials in total (' num2str(ntrialall-ntrial) ' had interval >1s, ' num2str(ntrialrej) ' had artifacts)'])
  disp(['accuracy fear   pred = ' num2str(fearpredacc*100,'%.1f') '%'])
  disp(['accuracy fear unpred = ' num2str(fearunpredacc*100,'%.1f') '%'])
  disp(['accuracy neut   pred = ' num2str(neutpredacc*100,'%.1f') '%'])
  disp(['accuracy neut unpred = ' num2str(neutunpredacc*100,'%.1f') '%'])
  disp(['accuracy        pred = ' num2str(predacc*100,'%.1f') '%'])
  disp(['accuracy      unpred = ' num2str(unpredacc*100,'%.1f') '%'])
  disp(['dprime          pred = ' num2str(norminv(fearpredacc)-norminv(1-neutpredacc),'%.2f')])
  disp(['dprime        unpred = ' num2str(norminv(fearunpredacc)-norminv(1-neutunpredacc),'%.2f')])
  disp(['accuracy        fear = ' num2str(fearacc*100,'%.1f') '%'])
  disp(['accuracy        neut = ' num2str(neutacc*100,'%.2f') ''])
  disp(['dprime               = ' num2str(norminv(fearacc)-norminv(1-neutacc),'%.2f')])
  disp(['rt mean(sd) fear   pred = ' num2str(fearpredrtm*1000,'%.0f') 'ms (' num2str(fearpredrts*1000,'%.0f') 'ms)'])
  disp(['rt mean(sd) fear unpred = ' num2str(fearunpredrtm*1000,'%.0f') 'ms (' num2str(fearunpredrts*1000,'%.0f') 'ms)'])
  disp(['rt mean(sd) neut   pred = ' num2str(neutpredrtm*1000,'%.0f') 'ms (' num2str(neutpredrts*1000,'%.0f') 'ms)'])
  disp(['rt mean(sd) neut unpred = ' num2str(neutunpredrtm*1000,'%.0f') 'ms (' num2str(neutunpredrts*1000,'%.0f') 'ms)'])
  disp(['rt mean(sd)        pred = ' num2str(predrtm*1000,'%.0f') 'ms (' num2str(predrts*1000,'%.0f') 'ms)'])
  disp(['rt mean(sd)      unpred = ' num2str(unpredrtm*1000,'%.0f') 'ms (' num2str(unpredrts*1000,'%.0f') 'ms)'])
  disp(['rt mean(sd)        fear = ' num2str(fearrtm*1000,'%.0f') 'ms (' num2str(fearrts*1000,'%.0f') 'ms)'])
  disp(['rt mean(sd)        neut = ' num2str(neutrtm*1000,'%.0f') 'ms (' num2str(neutrts*1000,'%.0f') 'ms)'])
  disp(['ntrial fear   pred = ' num2str(fearpredn)])
  disp(['ntrial fear unpred = ' num2str(fearunpredn)])
  disp(['ntrial neut   pred = ' num2str(neutpredn)])
  disp(['ntrial neut unpred = ' num2str(neutunpredn)])
  disp(' ')
  
  
end










%%%%%
% Fetch channel designations

% get details
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';

% set n's, loop, plot locations
nsubj = numel(info.subj);
chantype = [];
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(' ')
  disp(['patient ' currsubj])
  if isempty(info.sessionevents.(currsubj))
    disp('data not yet present')
    continue
  end
  
  % use only first session
  hdr = ft_read_header([info.datapath currsubj '/' info.sessiondata.(currsubj){1}]);
  label = setdiff(rmr_predfaceval_edftobesachannelnames(hdr.label),info.(currsubj).badchan);
  % remove numbers
  for ichan = 1:numel(label)
    label{ichan} = label{ichan}(isletter(label{ichan}));
  end
  chantype{isubj} = unique(label);
end
























%%%%%%%% old
% get details
info = rmr_predfaceval_info;

info.datapath = '/Users/roemer/Work/Data/tmpdata/pred_face_val/';

% adding back some thrown away patients
info.subj{end+1} = 'IR37';
info.sessionevents.IR37  = {'fp_task_13-Mar-2016_1','fp_task_13-Mar-2016_2','fp_task_13-Mar-2016_3','fp_task_13-Mar-2016_4'};

% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  disp(' ')
  disp(['patient ' currsubj])
  if isempty(info.sessionevents.(currsubj))
    disp('data not yet present')
    continue
  end
  
  
  % gather behavior
  predunpred  = [];
  fearneut  = [];
  percval  = [];
  accuracy = [];
  rt       = [];
  for    isess = 1   :numel(info.sessionevents.(currsubj))
    % set
    currsess = info.sessionevents.(currsubj){isess};
    
    % set fn
    fn = [info.datapath currsubj '/' currsess '.mat'];
    load(fn)
    predunpred  = [predunpred dat.cueType];
    fearneut    = [fearneut dat.faceValence];
    percval     = [percval dat.percValence];
    accuracy    = [accuracy dat.accuracy];
    rt          = [rt dat.rt];
    
    % cueType:       1 = predictive, 0 = nonpredictive
    % face valence:  1 = fear, 0 = neutral
    % percValence:   1 = fear, 0 = neural
    % accuracy:      1 = correct, 0 = incorrect
  end
  % based on rts
  rtmed = median(rt);
  rtmad = mad(rt,1);
  remind = (rt < (rtmed-6*rtmad)) | (rt > (rtmed+6*rtmad));
  rt(remind)         = [];
  accuracy(remind)   = [];
  fearneut(remind)   = [];
  predunpred(remind) = [];
  percval(remind)    = [];
  
  % compute/select
  ind = predunpred==1 & fearneut == 1;
  fearpredacc   = sum(accuracy(ind)) ./ sum(ind);
  fearpredrtm    = sum(rt(ind)) ./ sum(ind);
  %
  ind = predunpred==0 & fearneut == 1;
  fearunpredacc = sum(accuracy(ind)) ./ sum(ind);
  fearunpredrtm  = sum(rt(ind)) ./ sum(ind);
  %
  ind = predunpred==1 & fearneut == 0;
  neutpredacc   = sum(accuracy(ind)) ./ sum(ind);
  neutpredrtm    = sum(rt(ind)) ./ sum(ind);
  %
  ind = predunpred==0 & fearneut == 0;
  neutunpredacc = sum(accuracy(ind)) ./ sum(ind);
  neutunpredrtm  = sum(rt(ind)) ./ sum(ind);
  %
  ind = predunpred==1;
  predacc       = sum(accuracy(ind)) ./ sum(ind);
  predrtm        = sum(rt(ind)) ./ sum(ind);
  %
  ind = predunpred==0;
  unpredacc     = sum(accuracy(ind)) ./ sum(ind);
  unpredrtm      = sum(rt(ind)) ./ sum(ind);
  %
  ind = fearneut==1;
  fearacc       = sum(accuracy(ind)) ./ sum(ind);
  fearrtm        = sum(rt(ind)) ./ sum(ind);
  %
  ind = fearneut==0;
  neutacc       = sum(accuracy(ind)) ./ sum(ind);
  neutrtm        = sum(rt(ind)) ./ sum(ind);
  
  
  % spit it out
  disp([num2str(numel(rt)) ' trials in total'])
  disp(['accuracy fear   pred = ' num2str(fearpredacc*100,'%.1f') '%'])
  disp(['accuracy fear unpred = ' num2str(fearunpredacc*100,'%.1f') '%'])
  disp(['accuracy neut   pred = ' num2str(neutpredacc*100,'%.1f') '%'])
  disp(['accuracy neut unpred = ' num2str(neutunpredacc*100,'%.1f') '%'])
  disp(['accuracy        pred = ' num2str(predacc*100,'%.1f') '%'])
  disp(['accuracy      unpred = ' num2str(unpredacc*100,'%.1f') '%'])
  disp(['accuracy        fear = ' num2str(fearacc*100,'%.1f') '%'])
  disp(['accuracy        neut = ' num2str(neutacc*100,'%.1f') '%'])
  disp(['rt fear   pred = ' num2str(fearpredrtm*1000,'%.0f') 'ms'])
  disp(['rt fear unpred = ' num2str(fearunpredrtm*1000,'%.0f') 'ms'])
  disp(['rt neut   pred = ' num2str(neutpredrtm*1000,'%.0f') 'ms'])
  disp(['rt neut unpred = ' num2str(neutunpredrtm*1000,'%.0f') 'ms'])
  disp(['rt        pred = ' num2str(predrtm*1000,'%.0f') 'ms'])
  disp(['rt      unpred = ' num2str(unpredrtm*1000,'%.0f') 'ms'])
  disp(['rt        fear = ' num2str(fearrtm*1000,'%.0f') 'ms'])
  disp(['rt        neut = ' num2str(neutrtm*1000,'%.0f') 'ms'])
  
  
end





















% 
% 
% IR34: 243 usable trials in total (87 had interval >1s, 8 had artifacts)
% accuracy fear   pred = 69.9%
% accuracy fear unpred = 61.9%
% accuracy neut   pred = 83.1%
% accuracy neut unpred = 87.2%
% accuracy        pred = 76.0%
% accuracy      unpred = 75.3%
% dprime          pred = 1.48
% dprime        unpred = 1.44
% accuracy        fear = 67.2%
% accuracy        neut = 84.75
% dprime               = 1.47
% rt mean(sd) fear   pred = 1046ms (255ms)
% rt mean(sd) fear unpred = 983ms (241ms)
% rt mean(sd) neut   pred = 990ms (262ms)
% rt mean(sd) neut unpred = 954ms (270ms)
% rt mean(sd)        pred = 1020ms (259ms)
% rt mean(sd)      unpred = 968ms (256ms)
% rt mean(sd)        fear = 1025ms (251ms)
% rt mean(sd)        neut = 975ms (264ms)
% ntrial fear   pred = 83
% ntrial fear unpred = 42
% ntrial neut   pred = 71
% ntrial neut unpred = 47
%  
% 
% 
% IR35: 162 usable trials in total (38 had interval >1s, 98 had artifacts)
% accuracy fear   pred = 90.2%
% accuracy fear unpred = 93.8%
% accuracy neut   pred = 85.0%
% accuracy neut unpred = 79.3%
% accuracy        pred = 87.1%
% accuracy      unpred = 86.9%
% dprime          pred = 2.33
% dprime        unpred = 2.35
% accuracy        fear = 91.8%
% accuracy        neut = 83.15
% dprime               = 2.35
% rt mean(sd) fear   pred = 649ms (127ms)
% rt mean(sd) fear unpred = 630ms (114ms)
% rt mean(sd) neut   pred = 698ms (179ms)
% rt mean(sd) neut unpred = 768ms (179ms)
% rt mean(sd)        pred = 678ms (161ms)
% rt mean(sd)      unpred = 696ms (163ms)
% rt mean(sd)        fear = 641ms (121ms)
% rt mean(sd)        neut = 720ms (181ms)
% ntrial fear   pred = 41
% ntrial fear unpred = 32
% ntrial neut   pred = 60
% ntrial neut unpred = 29
% 
% 
% IR38: 168 usable trials in total (59 had interval >1s, 13 had artifacts)
% accuracy fear   pred = 97.9%
% accuracy fear unpred = 96.8%
% accuracy neut   pred = 98.2%
% accuracy neut unpred = 97.0%
% accuracy        pred = 98.1%
% accuracy      unpred = 96.9%
% dprime          pred = 4.14
% dprime        unpred = 3.72
% accuracy        fear = 97.5%
% accuracy        neut = 97.75
% dprime               = 3.96
% rt mean(sd) fear   pred = 760ms (143ms)
% rt mean(sd) fear unpred = 726ms (148ms)
% rt mean(sd) neut   pred = 723ms (88ms)
% rt mean(sd) neut unpred = 722ms (86ms)
% rt mean(sd)        pred = 740ms (117ms)
% rt mean(sd)      unpred = 724ms (119ms)
% rt mean(sd)        fear = 747ms (145ms)
% rt mean(sd)        neut = 723ms (87ms)
% ntrial fear   pred = 48
% ntrial fear unpred = 31
% ntrial neut   pred = 56
% ntrial neut unpred = 33
%  
% 
% 
% IR39: 226 usable trials in total (65 had interval >1s, 69 had artifacts)
% accuracy fear   pred = 88.2%
% accuracy fear unpred = 89.6%
% accuracy neut   pred = 87.5%
% accuracy neut unpred = 92.1%
% accuracy        pred = 87.9%
% accuracy      unpred = 90.7%
% dprime          pred = 2.34
% dprime        unpred = 2.67
% accuracy        fear = 88.8%
% accuracy        neut = 89.09
% dprime               = 2.45
% rt mean(sd) fear   pred = 691ms (170ms)
% rt mean(sd) fear unpred = 676ms (159ms)
% rt mean(sd) neut   pred = 663ms (112ms)
% rt mean(sd) neut unpred = 694ms (110ms)
% rt mean(sd)        pred = 677ms (143ms)
% rt mean(sd)      unpred = 684ms (139ms)
% rt mean(sd)        fear = 685ms (165ms)
% rt mean(sd)        neut = 674ms (112ms)
% ntrial fear   pred = 68
% ntrial fear unpred = 48
% ntrial neut   pred = 72
% ntrial neut unpred = 38
% 
% 
% 
% 
% 
% 
% 
% 
% IR41: 122 usable trials in total (39 had interval >1s, 182 had artifacts)
% accuracy fear   pred = 100.0%
% accuracy fear unpred = 100.0%
% accuracy neut   pred = 94.7%
% accuracy neut unpred = 95.0%
% accuracy        pred = 97.3%
% accuracy      unpred = 98.0%
% dprime          pred = Inf
% dprime        unpred = Inf
% accuracy        fear = 100.0%
% accuracy        neut = 94.83
% dprime               = Inf
% rt mean(sd) fear   pred = 1044ms (377ms)
% rt mean(sd) fear unpred = 1040ms (323ms)
% rt mean(sd) neut   pred = 892ms (344ms)
% rt mean(sd) neut unpred = 856ms (307ms)
% rt mean(sd)        pred = 965ms (366ms)
% rt mean(sd)      unpred = 965ms (326ms)
% rt mean(sd)        fear = 1042ms (351ms)
% rt mean(sd)        neut = 880ms (330ms)
% ntrial fear   pred = 35
% ntrial fear unpred = 29
% ntrial neut   pred = 38
% ntrial neut unpred = 20
% 
% 
% 
% IR43: 56 usable trials in total (16 had interval >1s, 134 had artifacts)
% accuracy fear   pred = 71.4%
% accuracy fear unpred = 70.0%
% accuracy neut   pred = 50.0%
% accuracy neut unpred = 30.0%
% accuracy        pred = 58.3%
% accuracy      unpred = 50.0%
% dprime          pred = 0.57
% dprime        unpred = 0.00
% accuracy        fear = 70.8%
% accuracy        neut = 43.75
% dprime               = 0.39
% rt mean(sd) fear   pred = 2983ms (1434ms)
% rt mean(sd) fear unpred = 2267ms (1642ms)
% rt mean(sd) neut   pred = 2463ms (1548ms)
% rt mean(sd) neut unpred = 3463ms (1681ms)
% rt mean(sd)        pred = 2665ms (1506ms)
% rt mean(sd)      unpred = 2865ms (1730ms)
% rt mean(sd)        fear = 2684ms (1532ms)
% rt mean(sd)        neut = 2776ms (1633ms)
% ntrial fear   pred = 14
% ntrial fear unpred = 10
% ntrial neut   pred = 22
% ntrial neut unpred = 10
% 
% 
% 
% IR44: 248 usable trials in total (92 had interval >1s, 19 had artifacts)
% accuracy fear   pred = 51.9%
% accuracy fear unpred = 64.4%
% accuracy neut   pred = 88.2%
% accuracy neut unpred = 91.3%
% accuracy        pred = 69.4%
% accuracy      unpred = 78.0%
% dprime          pred = 1.23
% dprime        unpred = 1.73
% accuracy        fear = 56.3%
% accuracy        neut = 89.34
% dprime               = 1.40
% rt mean(sd) fear   pred = 1437ms (549ms)
% rt mean(sd) fear unpred = 1458ms (506ms)
% rt mean(sd) neut   pred = 1357ms (507ms)
% rt mean(sd) neut unpred = 1289ms (548ms)
% rt mean(sd)        pred = 1398ms (529ms)
% rt mean(sd)      unpred = 1372ms (532ms)
% rt mean(sd)        fear = 1444ms (532ms)
% rt mean(sd)        neut = 1331ms (522ms)
% ntrial fear   pred = 81
% ntrial fear unpred = 45
% ntrial neut   pred = 76
% ntrial neut unpred = 46
% 
% 
























% IR34: 243 usable trials 
% accuracy        pred = 76.0%
% accuracy      unpred = 75.3%
% dprime          pred = 1.48
% dprime        unpred = 1.44
% rt mean(sd)        pred = 1020ms (259ms)
% rt mean(sd)      unpred = 968ms (256ms) 

% IR35: 162 usable trials 
% accuracy        pred = 87.1%
% accuracy      unpred = 86.9%
% dprime          pred = 2.33
% dprime        unpred = 2.35
% rt mean(sd)        pred = 678ms (161ms)
% rt mean(sd)      unpred = 696ms (163ms)

% IR38: 168 usable trials 
% accuracy        pred = 98.1%
% accuracy      unpred = 96.9%
% dprime          pred = 4.14
% dprime        unpred = 3.72
% rt mean(sd)        pred = 740ms (117ms)
% rt mean(sd)      unpred = 724ms (119ms)
%  
% 
% 
% IR39: 226 usable trials 
% accuracy        pred = 87.9%
% accuracy      unpred = 90.7%
% dprime          pred = 2.34
% dprime        unpred = 2.67
% rt mean(sd)        pred = 677ms (143ms)
% rt mean(sd)      unpred = 684ms (139ms)

% IR41: 122 usable trials 
% accuracy        pred = 97.3%
% accuracy      unpred = 98.0%
% dprime          pred = Inf
% dprime        unpred = Inf
% rt mean(sd)        pred = 965ms (366ms)
% rt mean(sd)      unpred = 965ms (326ms)

% IR44: 248 usable trials 
% accuracy        pred = 69.4%
% accuracy      unpred = 78.0%
% dprime          pred = 1.23
% dprime        unpred = 1.73
% rt mean(sd)        pred = 1398ms (529ms)
% rt mean(sd)      unpred = 1372ms (532ms)


