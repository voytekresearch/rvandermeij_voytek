function rmr_abcaudecog_behav



% get details
info = rmr_abcaudecog_info;


% set n's, loop, plot locations
nsubj = numel(info.subj);
for    isubj = 1    :nsubj
  
  % set
  currsubj = info.subj{isubj};
  switch currsubj
    case {'GP14','GP15','GP22','GP28','GP35','ST1','ST6','ST8'}
      fsample = 3051.76;
    case {'JH2','JH13'}
      fsample = 1000;
    otherwise
      error('woopsie')
  end
  
  
  % gather behavior
  attlrb     = [];
  target     = [];
  accuracy   = [];
  stimlrb    = [];
  rt         = [];
  ntrial     = 0;
  ntrialall  = 0;
  ntrialrej  = 0;
  for    isess = 1   :numel(info.session.(currsubj))
    % define trials
    cfg = [];
    cfg.subj    = currsubj;
    cfg.isess   = isess;
    cfg.fsample = fsample;
    trl = rmr_abcaudecog_definetrialsfrompos(cfg);
    
    % reject 
%     arttrl = info.(currsubj).artfctdef.visual.(info.session.(currsubj){isess}); 
%     trl    = cfgout.trl;
    
    % check if its empty :(
    if isempty(trl)
      warning([currdata ' contains no trials after artifact rejection'])
      continue
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
    
    attlrb     = [attlrb;      trl(:,4)];
    target     = [target;      trl(:,6)];
    accuracy   = [accuracy;    trl(:,7)];
    stimlrb    = [stimlrb;     trl(:,9)];
    rt         = [rt;          trl(:,8)];
    ntrial     = ntrial + size(trl,1);
    %     ntrialrej  = ntrialrej + (size(trlall,1) - size(trl,1) - sum(remind));
    %     ntrialall  = ntrialall + numel(remind);
  end
  %   rtmed = median(rt);
  %   rtmad = mad(rt,1);
  %   remind = (rt < (rtmed-6*rtmad)) | (rt > (rtmed+6*rtmad));
  %   rt(remind)         = [];
  %   attlrb(remind)   = [];
  %   target(remind)   = [];
  %   accuracy(remind) = [];
  %   ntrial = ntrial - sum(remind);
  
  % compute/select
  bool         = target == 1;
  alltarg      = sum(accuracy(bool)==1) ./ sum(bool);
  alltargrtm   = mean(rt(bool & accuracy == 1));
  alltargrts   = std(rt(bool & accuracy == 1));
  alltargn     = sum(bool);
  %
  bool = accuracy(accuracy==3 | accuracy==4) == 3;
  fa   = sum(bool) ./ numel(bool);
  nfa  = sum(bool);
  ncr  = sum(~bool);
  % 
  nattlstiml = sum(attlrb==1 & stimlrb==1);
  nattlstimr = sum(attlrb==1 & stimlrb==2);
  nattlstimb = sum(attlrb==1 & stimlrb==3);
  nattrstiml = sum(attlrb==2 & stimlrb==1);
  nattrstimr = sum(attlrb==2 & stimlrb==2);
  nattrstimb = sum(attlrb==2 & stimlrb==3);
  nattbstiml = sum(attlrb==3 & stimlrb==1);
  nattbstimr = sum(attlrb==3 & stimlrb==2);
  nattbstimb = sum(attlrb==3 & stimlrb==3);

  
  %
  % d prime = norminv(HR) - norminv(FAR)
  % c prime =  -(norminv(HR) + norminv(FAR))/2
  
  % spit it out
  disp([currsubj ': ' num2str(ntrial) ' trials, ' upper(info.(currsubj).leftrightelec) ' side electrodes'])
  disp(['hits targets        = ' num2str(alltarg*100,'%.1f') '%'])
  disp(['false alarms        = ' num2str(fa*100,'%.1f') '%'])
  disp(['rt mean(sd) targ    = ' num2str(alltargrtm,'%.0f') 'ms (' num2str(alltargrts,'%.0f') 'ms)'])
  disp(['ntrial targets      = ' num2str(alltargn)])
  disp(['ntrial false alarm  = ' num2str(nfa)])
  disp(['ntrial corr reject  = ' num2str(ncr)])
  disp(['ntrial att  left  - left  stim = ' num2str(nattlstiml)])
  disp(['ntrial att  left  - right stim = ' num2str(nattlstimr)])
  disp(['ntrial att  left  - bin   stim = ' num2str(nattlstimb)])
  disp(['ntrial att right  - left  stim = ' num2str(nattrstiml)])
  disp(['ntrial att right  - right stim = ' num2str(nattrstimr)])
  disp(['ntrial att right  - bin   stim = ' num2str(nattrstimb)])
  disp(['ntrial att   bin  - left  stim = ' num2str(nattbstiml)])
  disp(['ntrial att   bin  - right stim = ' num2str(nattbstimr)])
  disp(['ntrial att   bin  - bin   stim = ' num2str(nattbstimb)])
  disp(' ')
  
  
end




























