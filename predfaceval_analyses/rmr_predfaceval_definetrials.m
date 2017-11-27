function [trl,event] = rmr_predfaceval_definetrials(cfg)

% The function is used to segment the raw recordings using FieldTrip.
% The output of this function is a trl matrix, which describes the times of each trial in samples,
% and contains additional information for classifying trials (condition, hit/miss, rt, etc).
% The trl should be used as input for ft_preprocessing(cfg) by passing it as cfg.trl.
% The second output argument, event, is a structure containing most event information, with
% events coming from the diode.
%
% For all trials, t=0 is the onset of the CUE. The prestimulus period (tied to CUE) and
% poststimulus period (tied to FACE) can be defined using the options below.
%
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
%
% When using the trl as input for ft_preprocessing, columns 1-3 will be present in the
% output as data.sampleinfo. Columns 4-end will be present in data.trialinfo.
%
% The following options are necessary:
%   cfg.eventfile = string, filename of the .mat file output of PTB, linked to the specified datafile
%   cfg.datafile  = string, filename of the EDF/BESA file containing the raw recordings
%   cfg.diodechan = string, name of the analog channel containing the photo diode measurements, as specified in the EDF/BESA header
%   cfg.prestim   = scalar, latency in seconds before CUE ONSET
%   cfg.poststim  = scalar, latency in seconds after  FACE ONSET
% The following option(s) are optional:
%   cfg.threshold = 1x2 vector, thresholds for the upgoing (1) and downgoing (2) flank of diode events. For each cue/face event,
%                               the beginning and end points are found closest to these values. These values should be chosen for an individual
%                               session.
%
% This functions syncs the events recorded by the photo diode in the EDF/BESA file with the output in the PTB file.
% Syncing is performed on a per trial basis using the onset of the cue
% From the diode in the EDF/BESA file are taken:
%    - each cue onset (and therefore the start of the trial and t=0)
%    - each face onset (and therefore the end of the trial)
%    - each cue duration (additional row in the trl, see below)
%    - each face+isi duration (additional row in the trl, see below)
% From the PTB mat-file output are taken:
%    - all trial labels (valence, predictive or not, etc)
%    - response cue onset
%    - reaction time/sample
%
% To check whether the timings of the EDF/BESA and PTB are synced correctly, two deviations are checked:
%   1) deviations of experiment timeline after EDF/BESA and PTB are synced using cue of first trial (for tolerance see below)
%   2) per trial deviations of cue-face period (for tolerance see below)
%
% Trials on which the cue or face were on the screen for more than the specified frames are deleted.
% This is 200ms for the cue, and 50ms for the face+isi
%
%

% set defaults
cfg.threshold = ft_getopt(cfg, 'threshold', []);

% check for necessary options
if ~isfield(cfg, 'eventfile') || ~isfield(cfg, 'datafile')
  error('need to specify both cfg.eventfile and cfg.datafile')
end
if ~isfield(cfg, 'diodechan')
  error('need to specifiy name of channel containing diode information, use ft_databrowser to identify this channel')
end
if ~isfield(cfg, 'prestim') || ~isfield(cfg, 'poststim')
  error('need to specify cfg.prestim and cfg.poststim')
end
if isfield(cfg,'debugflg')
  debugflg = istrue(cfg.debugflg);
else
  debugflag = false;
end
if ~isempty(cfg.threshold) && numel(cfg.threshold)~=2
  error('cfg.threshold should contain two values')
end


% disp the input
datafnprts  = strsplit(cfg.datafile,{'\','/'});
eventfnprts = strsplit(cfg.eventfile,{'\','/'});
disp(['defining trials in ' datafnprts{end} ' and ' eventfnprts{end}])

% read the header information
hdr = ft_read_header(cfg.datafile);


%%%%%%%%%%%
%%% Diode event detection
% read the events from the data
chanindx    = find(ismember(hdr.label, ft_channelselection(cfg.diodechan, hdr.label)));
% readin event channel to use later
diodesig = ft_read_data(cfg.datafile);
diodesig = diodesig(chanindx,:); % bypass reading error I can't fix right now
% detect up and down going flanks
event = ft_read_event(cfg.datafile, 'chanindx', chanindx, 'detectflank', 'both', 'threshold', '(3/2)*nanmedian');
% estimate DAC noise/minimal step in diode signal, to use for correction
diodediff = sort(abs(diff(diodesig)),'ascend');
diodediff(diodediff==0) = [];
diodestep = mean(diodediff(1:round(hdr.Fs*0.100))); % take 100ms of the minimal diffs, should be enough
% prune events, which is necessary due to DAC noise (the diode signal can flip one DAC value back and forth, which can cause faulty event detection)
falseeventind = [];
for ievent = 1:numel(event)
  switch event(ievent).type
    case [cfg.diodechan '_up']
      % when up, 75% of the next non-constant 20 samples should agree
      futsmp = diodesig(event(ievent).sample+1:(event(ievent).sample+20));
      nconst = sum(diff(futsmp)==0);
      if sum(sign(diff(futsmp)) == 1) < ((20-nconst)*.75)
        falseeventind = [falseeventind ievent];
      end
    case [cfg.diodechan '_down']
      % when down, 75% of the previous non-constant 50 samples should agree (downgoing flank flat)
      futsmp = diodesig((event(ievent).sample-50):(event(ievent).sample-1));
      nconst = sum(diff(futsmp)==0);
      if sum(sign(diff(futsmp)) == -1) <((50-nconst)*.75)
        falseeventind = [falseeventind ievent];
      end
  end
end
event(falseeventind) = [];

% further prune events by removing accidental duplicates: remove ups and downs from the same flank, keeping the first up and the first down
upind      = find(strcmp({event.type},[cfg.diodechan '_up']));
upinddelay = diff([NaN event(upind).sample]) ./ hdr.Fs;
upindrep   = [NaN diff(upind)] == 1 & upinddelay<(1/60); % replicate must occur within the time of a frame at 60Hz
event(upind(upindrep)) = [];
downind      = find(strcmp({event.type},[cfg.diodechan '_down']));
downinddelay = diff([NaN event(downind).sample]) ./ hdr.Fs;
downindrep   = [NaN diff(downind)] == 1 & downinddelay<(1/60); % replicate must occur within the time of a frame at 60Hz
event(downind(downindrep)) = [];
% report on number of events removed
disp([num2str(sum([numel(falseeventind) sum(upindrep) sum(downindrep)])) ' single diode events (' num2str(sum([sum(upindrep) sum(downindrep)])) ' replicates) were incorrectly detected and skipped'])

% now that we can assume only real up/down events are found, increase their accuracy
for ievent = 1:numel(event)
  switch event(ievent).type
    case [cfg.diodechan '_up']
      if isempty(cfg.threshold)
        % search backward in time with a sliding window until the start of the up-flank has been found
        % (diode signal rise speed increases with time)
        found = 0;
        lastsmp = event(ievent).sample;
        while ~found
          nsmp = 40;
          smpwin = (lastsmp-nsmp-1):(lastsmp);
          winmean = mean(diodesig(smpwin));
          if (winmean-diodesig(lastsmp))<(-diodestep*2.5) % difference between window and last sample should be at least 2 steps (2.5 to account for small variotions in stepsize)
            % continue search
            lastsmp = lastsmp - 1;
          else
            % end found
            lastsmp = lastsmp + 1; % previous was correct one
            found = true;
            event(ievent).sample = lastsmp;
          end
        end
      else
        % find first sample after crossing threshold, searching backwards maximally a 50ms
        oldsmp     = event(ievent).sample;
        nsmpsearch = round(0.050*hdr.Fs);
        searchsig  = diodesig((oldsmp-nsmpsearch):oldsmp); 
        ind        = find(diff(sign(searchsig-cfg.threshold(1))),1,'last');
        newsmp     = oldsmp - (nsmpsearch - ind);
        event(ievent).sample = newsmp;
      end
      
    case [cfg.diodechan '_down']
      if isempty(cfg.threshold)
        % search backward in time with a sliding window until the start of the down-flank has been found
        % (diode signal decay speed decreases with time)
        found = 0;
        lastsmp = event(ievent).sample;
        while ~found
          nsmp = 40;
          smpwin = (lastsmp-nsmp-1):(lastsmp);
          winmean = mean(diodesig(smpwin));
          if (winmean-diodesig(lastsmp))>(diodestep*2.5) % difference between window and last sample should be at least 2 steps (2.5 to account for small variotions in stepsize)
            % continue search
            lastsmp = lastsmp - 1;
          else
            % end found
            lastsmp = lastsmp + 1; % previous was correct one
            found = true;
            event(ievent).sample = lastsmp;
          end
        end
      else
        % find first sample before crossing threshold, searching backwards maximally a 50ms
        oldsmp     = event(ievent).sample;
        nsmpsearch = round(0.050*hdr.Fs);
        searchsig  = diodesig((oldsmp-nsmpsearch):oldsmp); 
        ind        = find(diff(sign(searchsig-cfg.threshold(2))),1,'last');
        newsmp     = oldsmp - (nsmpsearch - ind);
        event(ievent).sample = newsmp;
      end
      
  end
end
% convert up and down events into events with a duration in seconds starting at the up event
preparseevent = event;
event = [];
count = 0;
for ievent = 1:numel(preparseevent)
  if strcmp(preparseevent(ievent).type,[cfg.diodechan '_up'])
    if numel(preparseevent)>=(ievent+1) && strcmp(preparseevent(ievent+1).type,[cfg.diodechan '_down'])
      count = count+1;
      begsample = preparseevent(ievent).sample;
      endsample = preparseevent(ievent+1).sample;
      event(count).sample   = begsample;
      event(count).duration = (endsample-begsample+1)./hdr.Fs;
      % determine event type using broad margins
      if (event(count).duration > (0.2-(2.5/60)) && event(count).duration <  (0.2+(2.5/60))) % defined length is 200ms with errors of +/- 2.5/60 (2 refresh ticks)
        event(count).type = 'cue';
      elseif (event(count).duration > (0.05-(2.5/60)) && event(count).duration < (0.05+(2.5/60)))
        event(count).type = 'face+isi';
      else
        event(count).type = 'other';
      end
    end
  end
end
disp([num2str(numel(event)*2) ' single diode events were paired'])
disp([num2str(numel(preparseevent)-(numel(event)*2)) ' single diode events were unpaired and skipped'])

% find indices of full trial events and remove others
precleanevent = event;
% first remove other events
othereventind = find(strcmp({event.type},'other'));
event(othereventind) = [];
disp([num2str(numel(othereventind)) ' events were not cue/face+isi and removed'])
% then, only keep events where cue is followed by a face+isi
keepind = [];
for ievent = 1:(numel(event)-rem(numel(event),2))
  if strcmp(event(ievent).type,'cue')
    if strcmp(event(ievent+1).type,'face+isi')
      keepind = [keepind ievent ievent+1];
    end
  end
end
disp([num2str(numel(event)-numel(keepind)) ' cue/face+isi events occured in isolation and were removed'])
event  = event(keepind);
ntrial = numel(event)/2;

%%%%%
% diode debug plotting
% plot photo diode signal with detected event onsets/offsets
if debugflg
  figure('numbertitle','off','name',[datafnprts{end} ' and ' eventfnprts{end}])
  subplot(2,1,1)
  hold on
  plot((1:numel(diodesig))./hdr.Fs,diodesig)
  plot([preparseevent.sample]./hdr.Fs,diodesig([preparseevent.sample]),'sc','markerfacecolor','cyan')
  plot([event.sample]./hdr.Fs,diodesig([event.sample]),'marker','v','color','r')
  if isempty(othereventind) % happens often
    legend('diode signal','detected up/down events','final cue/face+isi events')
  else
    plot([precleanevent(othereventind).sample]./hdr.Fs,diodesig([precleanevent(othereventind).sample]),'vgr')
    legend({'diode signal','detected up/down events','final cue/face+isi events','removed cue/face+isi events'})
  end
  if ~isempty(cfg.threshold)
    x = [1 numel(diodesig)]./hdr.Fs;
    y = [cfg.threshold(1) cfg.threshold(1)];
    line(x,y,'color',[.5 .5 .5],'linestyle','--')
    y = [cfg.threshold(2) cfg.threshold(2)];
    line(x,y,'color',[.5 .5 .5],'linestyle','--')
  end
  xlabel('time(s)')
  ylabel('diode signal strength')
  title(['diode events of ' datafnprts{end}],'interpreter','none')
end
% diode debug plotting
%%%%%

% act upon found events
if ~isempty(event)
  cueind  = find(strcmp({event.type},'cue'));
  faceind = find(strcmp({event.type},'face+isi'));
  cuefaceoffs = (([event(faceind).sample]-[event(cueind).sample]) ./ hdr.Fs) - [event(cueind).duration];
  cuedur  = [event(cueind).duration];
  facedur = [event(faceind).duration];
  disp(['found ' num2str(ntrial) ' trials using EDF/BESA diode recording'])
  disp(['EDF/BESA diode        cue duration mean(sd) = ' num2str(mean(cuedur))      's (' num2str(std(cuedur))      's)  min = ' num2str(min(cuedur))      's max = ' num2str(max(cuedur)) 's'])
  disp(['EDF/BESA diode   face+isi duration mean(sd) = ' num2str(mean(facedur))     's (' num2str(std(facedur))     's)  min = ' num2str(min(facedur))     's max = ' num2str(max(facedur)) 's'])
  disp(['EDF/BESA diode cue-face+isi offset mean(sd) = ' num2str(mean(cuefaceoffs)) 's (' num2str(std(cuefaceoffs)) 's)  min = ' num2str(min(cuefaceoffs)) 's max = ' num2str(max(cuefaceoffs)) 's'])
else
  error('no trials were found in EDF/BESA diode recording')
end
%%%%%%%%%%%


%%%%%%%%%%%
% obtain events from PTB mat-file
% read output from PTB mat-file
ptb = load(cfg.eventfile);
% check whether number of detected trials in PTB and EDF/BESA correspond
ptbntrial = numel(ptb.dat.cueType);
disp(['found ' num2str(ptbntrial) ' trials in PTB output'])
% deal with number of trials in PTB mat-file being different from that in EDF/BESA file
if ntrial~=ptbntrial
  if ntrial<ptbntrial
    % two possible causes:
    % 1) too few events detect, event diode detection failed
    % 2) recording was switched on too late, or switched off too early (*sigh*)
    % Rather than allowing an exception to occur (in the case of cause 2) for every subject, I'm only gonna allow it for those
    % in which I know it occurs, which is safer
    
    % parse file name to identify dataset
    [path sessname ext] = fileparts(cfg.datafile);
    switch sessname
      case '2016051512_0005' % a session from IR41 where the recording was turned on too late...
        % the first trial was not recorded in the diode, cut it out of all dat fields that are used...
        ptb.dat.cueType      = ptb.dat.cueType(2:end);
        ptb.dat.faceValence  = ptb.dat.faceValence(2:end);
        ptb.dat.percValence  = ptb.dat.percValence(2:end);
        ptb.dat.accuracy     = ptb.dat.accuracy(2:end);
        ptb.dat.rt           = ptb.dat.rt(2:end);
        ptb.dat.cueTargetInterval  = ptb.dat.cueTargetInterval(2:end);
        ptb.dat.interTrialInterval = ptb.dat.interTrialInterval(2:end);
        ptb.dat.faceIdentity = ptb.dat.faceIdentity(2:end);
        if isfield(ptb.dat,'frames') % time stamps are present
          remind = 1:find(diff(ptb.dat.frameID==10)==-1,1); % removes trial defined by last response frame of first trial
          ptb.dat.frameID(remind) = [];
          ptb.dat.frames(remind) = [];
        end
        ptbntrial = 59;
        disp('exception for IR42: 2016051512_0005 - recording switched on too late, using last 59 trials from PTB output')
      otherwise
        error('different number of trials detected in EDF/BESA file and PTB output, evaluate diode event detection or add exception')
    end
  else
    % too many events detected, diode event detection failed
    error('different number of trials detected in EDF/BESA file and PTB output')
  end
end
% change some of the coding
ptb.dat.cueType = ~ptb.dat.cueType + 1;          % is  1 = predictive, 0 = nonpredictive, change to 1 = predictive, 2 = unpredictive
ptb.dat.faceValence = ~ptb.dat.faceValence + 1;  % is  1 = fear, 0 = neutral,             change to 1 = fear, 2 = neutral
ptb.dat.percValence = ~ptb.dat.percValence +1;   % is  1 = fear, 0 = neural,              change to 1 = fear, 2 = neutral


% extract information from PTB mat-file to use for checking syncing accuracy
if isfield(ptb.dat,'frames') % time stamps are present
  hastimestamps = true;
  % fetch cue onset from timestamps
  cueframes      = ptb.dat.frameID==2; % 2 = cue
  ptbcueonset    = ptb.dat.frames(diff(cueframes)==1);
  % fetch face onset from timestamps
  faceframes     = ptb.dat.frameID==4; % 4 = target
  ptbfaceonset   = ptb.dat.frames(diff(faceframes)==1);
  % fetch response cue onset from timestamps
  faceframes     = ptb.dat.frameID==10; % 10 = response screen
  ptb.respcueons = ptb.dat.frames(diff(faceframes)==1);
  ptb.respcueons = ptb.respcueons - ptbcueonset; % these should be relative from cue onset (t=0)
  
else % time stamps are missing... use workaround to get inaccurate representation of timeline
  hastimestamps = false;
  % reconstruct ptbcueonset, used for syncing
  ifi = 1/60; % assume 60Hz refresh rate
  bs  = ones(1,ntrial);
  timeline = [ptb.dat.interTrialInterval;...   % iti
    round(bs * .2 ./ ifi) * ifi;...            % cue
    ptb.dat.cueTargetInterval;...              % delay
    ptb.dat.rt;...                             % rt (rt is counted from face onset onwards, so includes face+soa+4*mask
    round(bs .* 0.2 ./ ifi) * ifi];            % delay after response
  timeline = timeline(:);
  timeline = cumsum(timeline);
  % get cue onset
  ptbcueonset = timeline(2:5:end)';
  % start at zero
  ptbcueonset = ptbcueonset - ptbcueonset(1);
  % correct for some of the cumulative timing errors by stretching ptbcueonset (these errors are due to dropped frames)
  % get optimal strecht factor using regression
  reccueonset = [event(1:2:end).sample] ./ hdr.Fs;
  reccueonset = reccueonset - reccueonset(1); % start timeline at 0
  scalefac = reccueonset / ptbcueonset;
  ptbcueonset  = ptbcueonset .* scalefac;
  
  % construct faceonset from cueonset
  ptbfaceonset = ptbcueonset + ptb.dat.cueTargetInterval + (round(.2 ./ ifi) * ifi);
  
  % construct respcueons
  ptb.respcueons = NaN(1,ntrial);
end
%%%%%%%%%%%



%%%%%%%%%%%
%%% Perform syncing and timing check
ptbfaceonset = ptbfaceonset - ptbcueonset(1); % start timeline at 0
ptbcueonset  = ptbcueonset - ptbcueonset(1); % start timeline at 0
reccueonset  = [event(1:2:end).sample] ./ hdr.Fs;
recfaceonset = [event(2:2:end).sample] ./ hdr.Fs;
recfaceonset = recfaceonset - reccueonset(1); % start timeline at 0
reccueonset  = reccueonset - reccueonset(1); % start timeline at 0
% syncing check 1 - calc syncing error after aligning to cue of first trial
cosyncerror = ptbcueonset-reccueonset;
cosyncerror = cosyncerror * 1000;
% syncing check 2 - calc syncing error between cue and face onsets
ptbcfonsetdiff = ptbfaceonset - ptbcueonset;
reccfonsetdiff = ([event(2:2:end).sample]-[event(1:2:end).sample]) ./ hdr.Fs;
cfodsyncerror = ptbcfonsetdiff-reccfonsetdiff;
cfodsyncerror = cfodsyncerror*1000;

% timing check: remove trials where cue/face+isi durations were wrong
cuedur     = [event(1:2:end).duration];
facedur    = [event(2:2:end).duration];
cuefaceint = ([event(2:2:end).sample]-[event(1:2:end).sample]+1) ./ hdr.Fs - cuedur;
% use 5ms, to capture not only frame errors, but also serious diode errors (not, at 5khz, 5ms is an error of 25 samples)
remind = (cuedur < (.200-.005)) | (cuedur > (.200+.005)) | (facedur < (.050-.005)) | (facedur > (.050+.005)) | (abs(cuefaceint-ptb.dat.cueTargetInterval)>0.005);
% instead of removing the trail from all possible fields, remove them from the selection
trialind = 1:ptbntrial;
trialind(remind) = [];
ntrial   = numel(trialind);

%%%%%
% BESA/EDF - PTB syncing debug plotting
if debugflg
  subplot(2,1,2)
  hold on
  lincol = lines(ptbntrial);
  % plot diode cue/face+isi
  lincol(remind,:) = repmat([.8 .8 .8],[sum(remind) 1]); % grey out removed trials
  reccfpairs = [reccueonset; recfaceonset];
  for itrial = 1:ptbntrial
    plot(reccfpairs(:,itrial),[1.05 1.05],'marker','^','color',lincol(itrial,:))
  end
  % plot ptb cue/face+isi
  ptbcfpairs = [ptbcueonset; ptbfaceonset];
  for itrial = 1:ptbntrial
    plot(ptbcfpairs(:,itrial),[1 1],'marker','v','color',lincol(itrial,:))
  end
  % plot connecting lines
  for itrial = 1:ptbntrial
    plot([reccfpairs(1,itrial) ptbcfpairs(1,itrial)],[1.05 1],'color',lincol(itrial,:))
    plot([reccfpairs(2,itrial) ptbcfpairs(2,itrial)],[1.05 1],'color',lincol(itrial,:))
  end
  set(gca,'ylim',[0.8 1.2],'ytick',[1 1.05],'yticklabel',{'PTB','BESA/EDF'},'xlim',[-2 max([reccfpairs(:); ptbcfpairs(:)])+2])
  xlabel('time(s)')
  title('cue/face+isi timing per trial of BESA/EDF vs PTB (grey = removed due to timing errors)')
  legend('trial 1', 'trial 2', 'trial 3', 'etc')
  drawnow % flush
end
%%%%%

% act upon calculated syncing errors
% syncing check 1 - syncing error after aligning to cue of first trial
disp(['experiment-wide recording-ptb timing offset: max = ' num2str(max(abs(cosyncerror))) 'ms  med(sd) = ' num2str(median(cosyncerror)) 'ms (' num2str(std(cosyncerror)) 'ms)'])
if hastimestamps
  if max(abs(cosyncerror))>200
    error('severe error (>200ms) detected in syncing recording and PTB experiment-wide time-axis')
  end
else
  if     max(abs(cosyncerror))>200 && max(abs(cosyncerror))<2000
    warning('error (>200ms & <2000ms) detected in syncing recording and PTB experiment-wide time-axis (note, no time-stamps were found)')
  elseif max(abs(cosyncerror))>2000
    error('severe error (>2000ms) detected in syncing recording and PTB experiment-wide time-axis (note, no time-stamps were found)')
  end
end
% syncing check 2 - syncing error between cue and face onsets
disp(['trial-specific recording-ptb timing offset of cue-face delay: max = ' num2str(max(abs(cfodsyncerror))) 'ms  med(sd) = ' num2str(median(cfodsyncerror)) 'ms (' num2str(std(cfodsyncerror)) 'ms)'])
if hastimestamps
  if max(abs(cfodsyncerror))>20
    error('severe error (>20ms) detected in syncing recording and PTB trial-by-trial time-axis') % 1 refresh tick is 1/60 = 16.7ms
  end
else
  if      max(abs(cfodsyncerror))>20 && max(abs(cfodsyncerror))<100
    warning('error (>20ms & <100ms) detected in syncing recording and PTB trial-by-trial time-axis (note, no time-stamps were found)') % 1 refresh tick is 1/60 = 16.7ms
  elseif  max(abs(cfodsyncerror))>100
    error('severe error (>100ms) detected in syncing recording and PTB trial-by-trial time-axis (note, no time-stamps were found)') % 1 refresh tick is 1/60 = 16.7ms
  end
end
% sync was succesful
disp('syncing deviations are within tolerance')
% timing check
disp([num2str(ptbntrial-ntrial) ' trials had cue/face+isi duration timing errors that exceeded 5ms and were removed'])
if (ptbntrial-ntrial)>5
  warning('more than 5 trials were removed due cue/face+isi duration timing errors that exceeded 5ms')
end
%%%%%%%%%%%



% syncing is accurate enough, proceed with building trial definition
disp(['creating trl matrix containing ' num2str(ntrial) ' trials with mean(SD) cue->face onset differences of ' num2str(mean(reccfonsetdiff)) 's (' num2str(std(reccfonsetdiff)) 's)'])

% see documentation above for individual columns
trl = [];
for itrial = trialind
  curreventind = [(itrial*2)-1 (itrial*2)];
  
  % get trl samples
  begsample = event(curreventind(1)).sample -round(cfg.prestim  * hdr.Fs); % cue onset  - prestim
  endsample = event(curreventind(2)).sample +round(cfg.poststim * hdr.Fs); % face onset + poststim
  offset    = -round(cfg.prestim*hdr.Fs);
  % get trial info from PTB
  predunpred = ptb.dat.cueType(itrial);      % 1 = predictive, 2 = unpredictive
  fearneut   = ptb.dat.faceValence(itrial);  % 1 = fear, 2 = neutral
  respval    = ptb.dat.percValence(itrial);  % 1 = fear, 2 = neutral
  hitmiss    = fearneut == respval;          % 1 = hit, 0 = miss
  rt         = ptb.dat.rt(itrial);           % in seconds
  faceindex  = ptb.dat.faceIdentity(itrial);
  
  % get various latencies
  faceonset  = ((event(curreventind(2)).sample)-(event(curreventind(1)).sample)+1) ./ hdr.Fs; % t=0 is sample=1(-offset)
  respcueons = ptb.respcueons(itrial);
  % respcueons could be on a timepoint that doesn't match a sample, find closest one
  respcueons = round(respcueons .* hdr.Fs) ./ hdr.Fs;
  
  % grab duration of cue and face+isi
  cuediodur  = event(curreventind(1)).duration;
  facediodur = event(curreventind(2)).duration;
  
  % put all in trl
  trl(end+1,:) = [begsample endsample offset predunpred fearneut hitmiss rt faceonset respcueons faceindex cuediodur facediodur];
end
% remove trials that fall outside of bounds
if any(trl(:,1)<1)
  warning([num2str(sum(trl(:,1)<1)) ' trial(s) start before the recording due to cfg.prestim and were removed' ])
  trl(trl(:,1)<1,:) = [];
end
if any(trl(:,2)>hdr.nSamples)
  warning([num2str(sum(trl(:,2)>hdr.nSamples)) ' trial(s) end after the recording due to cfg.poststim and were removed' ])
  trl(trl(:,2)>hdr.nSamples,:) = [];
end



% modify event structure for output
for ievent = 1:numel(event)
  % add whether trials were removed due to timing errors
  if any(ievent == (trialind*2-1)) || any(ievent == (trialind*2))
    event(ievent).timingerror = false;
  else
    event(ievent).timingerror = true;
  end
  
  % add trial information from PTB
  if mod(ievent,2)==0% even
    currtrial = ievent/2;
  else % odd
    currtrial = ceil(ievent/2);
  end
  event(ievent).predunpred = ptb.dat.cueType(currtrial);      % 1 = predictive, 2 = unpredictive
  event(ievent).fearneut   = ptb.dat.faceValence(currtrial);  % 1 = fear, 2 = neutral
  event(ievent).respval    = ptb.dat.percValence(currtrial);  % 1 = fear, 2 = neutral
  event(ievent).hitmiss    = event(ievent).fearneut == event(ievent).respval;  % 1 = hit, 0 = miss
  event(ievent).rt         = ptb.dat.rt(currtrial);           % in seconds
  event(ievent).faceindex  = ptb.dat.faceIdentity(currtrial);
end
% set event.type to value and set type to duide per FT convention
for ievent = 1:numel(event)
  event(ievent).value  = event(ievent).type;
  event(ievent).type  = 'diode';
end


























