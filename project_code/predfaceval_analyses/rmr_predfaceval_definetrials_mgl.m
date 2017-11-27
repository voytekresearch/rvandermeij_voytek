function trl = rmr_irvfacepred_definetrials(cfg)

% The function is used to segment the raw recordings using FieldTrip. 
% The output of this function is a trl matrix, which describes the times of each trial in samples, 
% and contains additional information for classifying trials (condition, hit/miss, rt, etc). 
% The trl should be used as input for ft_preprocessing(cfg) by passing it as cfg.trl. 
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
%  8) faceonset   - onset of face stimulus in seconds
%  9) respcueons  - onset of response cue in seconds 
% 10) faceindex   - index of shown face as defined in MGL task script
% 11) cuediodur   - duration of cue on the screen as measured by the photo diode
% 12) facediodur  - duration of face+mask on the screen as measured by the photo diode
% 
% When using the trl as input for ft_preprocessing, columns 1-3 will be present in the 
% output as data.sampleinfo. Columns 4-end will be present in data.trialinfo.
%
% The following options are necessary:
%   cfg.eventfile = string, filename of the .mat file output of MGL, linked to the specified datafile
%   cfg.datafile  = string, filename of the EDF file containing the raw recordings
%   cfg.diodechan = string, name of the analog channel containing the photo diode measurements, as specified in the EDF header
%   cfg.prestim   = scalar, latency in seconds before CUE ONSET 
%   cfg.poststim  = scalar, latency in seconds after  FACE ONSET 
%
% This functions syncs the events recorded by the photo diode in the EDF file with the output in the MGL file. 
% Syncing is performed on a per trial basis using the onset of the cue
% From the diode in the EDF file are taken:
%    - each cue onset (and therefore the start of the trial and t=0)
%    - each face onset (and therefore the end of the trial)
%    - each cue duration (additional row in the trl, see below)
%    - each face+mask duration (additional row in the trl, see below)
% From the MGL output is taken:
%    - all trial labels (valence, predictive or not, etc)
%    - response cue onset
%    - reaction time/sample
%
% To check whether the timings of the EDF/MGL are synced correctly, two deviations are checked:
%   1) deviations of experiment timeline after EDF/GML are synced using cue of first trial (tolerance is max(abs) = 100ms)
%   2) per trial deviations of cue-face period (tolerance is max(abs) = 20ms, 1 refresh tick is 1/60 = 16.7ms)
%
% 


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

% read the header information
hdr = ft_read_header(cfg.datafile);

% read the events from the data
chanindx      = find(ismember(hdr.label, ft_channelselection(cfg.diodechan, hdr.label)));
detectflank   = 'both'; % detect up and down flanks
threshold     = '(3/2)*nanmedian';
event         = ft_read_event(cfg.datafile, 'chanindx', chanindx, 'detectflank', detectflank, 'threshold', threshold);
% convert up and down events into events with a duration in seconds starting at the up event
oldevent = event;
event = [];
count = 0;
for ievent = 1:numel(oldevent)
  if strcmp(oldevent(ievent).type,[cfg.diodechan '_up'])
    if numel(oldevent)>=(ievent+1) && strcmp(oldevent(ievent+1).type,[cfg.diodechan '_down'])
      count = count+1;
      begsample = oldevent(ievent).sample;
      endsample = oldevent(ievent+1).sample;
      event(count).sample   = begsample;
      event(count).duration = (endsample-begsample+1)./hdr.Fs;
      % determine event type using broad margins
      if (event(count).duration > (0.2-.025) && event(count).duration <  (0.2+.025)) % defined length is 200ms with errors of +/- 1/60 (refresh ticks)
        event(count).type = 'cue';
      elseif (event(count).duration > (0.05-0.025) && event(count).duration < (0.05+0.025))
        event(count).type = 'face+mask';
      else
        event(count).type = 'other';
      end
    end
  end
end
if (numel(oldevent)-(numel(event)*2))~=0
  warning([num2str(numel(oldevent)-(numel(event)*2)) ' diode events were incorrectly detected and skipped'])
end

% %%% DEBUG PLOTTING
% % plot photo diode signal with detected event onsets/offsets
% dat = ft_read_data(cfg.datafile);
% dat = dat(chanindx,:); % bypass error I can't fix right now
% figure
% hold on
% plot(dat)
% plot([event.sample],dat([event.sample]),'marker','*')
% %%% DEBUG PLOTTING


% keep only events where a cue is followed by a face+mask 
keepind = [];
for ievent = 1:numel(event)
  if strcmp(event(ievent).type,'cue')
    if strcmp(event(ievent+1).type,'face+mask')
      keepind = [keepind ievent ievent+1];
    end
  end
end
event  = event(keepind);
ntrial = numel(event)/2;
if isempty(event) || numel(event)==1
  error('no cue events detected that are followed by a face+mask event: no trials found')
end
disp(['found ' num2str(ntrial) ' trials using EDF diode recording'])

% read events from MGL output using MGL function copied as subfunction below
load(cfg.eventfile) % contains myscreen/stimulus/task
mgltrialinfo = getTaskParameters(myscreen,task);
disp(['found ' num2str(mgltrialinfo.nTrials) ' trials using MGL output'])

% throw error if number of trials in MGL output is different from that in EDF file
if ntrial~=mgltrialinfo.nTrials
  error('different number of trials detected in EDF file and MGL output')
end

%%%%%
% mgltrialinfo.trials(itrial).segtime contains trial and timing information from mgl
% Each element of it refers to:
%  (1) iti
%  (2) cue
%  (3) delay
%  (4) target 
%  (5) soa
%  (6) mask1
%  (7) mask2
%  (8) mask3
%  (9) mask4
% (10) response
% (11) response delay
%%%%%

%%%%%%%%%%%
%%% Perform syncing check
% check 1 - syncing error after aligning to cue of first trial
mglcueonset  = zeros(1,ntrial);
for itrial = 1:ntrial
  mglcueonset(itrial) = mgltrialinfo.trials(itrial).segtime(2);
end
mglcueonset = mglcueonset - mglcueonset(1); % start timeline at 0
edfcueonset = [event(1:2:end).sample] ./ hdr.Fs;
edfcueonset = edfcueonset - edfcueonset(1); % start timeline at 0
syncerror = mglcueonset-edfcueonset;
syncerror = syncerror * 1000;
disp(['experiment-wide edf-mgl timing offset in milisecond: max = ' num2str(max(abs(syncerror))) '  med(sd) = ' num2str(median(syncerror)) ' (' num2str(std(syncerror)) ')'])
if max(abs(syncerror))>100
  error('severe error (>100ms) detected in syncing EDF and MGL experiment-wide time-axis')
end
% check 2 - syncing error between cue and face onsets
mglcfonsetdiff = zeros(1,ntrial);
for itrial = 1:ntrial
  mglcfonsetdiff(itrial) = mgltrialinfo.trials(itrial).segtime(4)-mgltrialinfo.trials(itrial).segtime(2);
end
edfcfonsetdiff = ([event(2:2:end).sample]-[event(1:2:end).sample]) ./ hdr.Fs;
syncerror = mglcfonsetdiff-edfcfonsetdiff;
syncerror = syncerror*1000;
disp(['trial-specific edf-mgl timing offset of cue-face delay in milisecond: max = ' num2str(max(abs(syncerror))) '  med(sd) = ' num2str(median(syncerror)) ' (' num2str(std(syncerror)) ')'])
if max(abs(syncerror))>20
  error('error (>20ms) detected in syncing EDF and MGL trial-by-trial time-axis') % 1 refresh tick is 1/60 = 16.7ms
end
%%%%%%%%%%%

% syncing is accurate enough, proceed with building trial definition
disp(['syncing deviations within tolerance, creating trl matrix containing ' num2str(ntrial) ' trials with mean(SD) cue->face onset differences of ' num2str(mean(edfcfonsetdiff)) 's (' num2str(std(edfcfonsetdiff)) ')'])

% see documentation above for individual columns
trl = [];
eventind = 1:2:(ntrial*2); % event is serialized structure of 2 events per trial
for itrial = 1:ntrial
  curreventind = eventind(itrial);
 
  % get trl samples
  begsample = event(curreventind).sample   -round(cfg.prestim  * hdr.Fs); % cue onset  - prestim
  endsample = event(curreventind+1).sample +round(cfg.poststim * hdr.Fs); % face onset + poststim
  offset    = -round(cfg.prestim*hdr.Fs);
  % get trial info from MGL
  predunpred = mgltrialinfo.randVars.cuetype(itrial);     % 1 = predictive, 2 = unpredictive
  fearneut   = mgltrialinfo.randVars.facevalence(itrial); % 1 = fear, 2 = neutral
  respval    = mgltrialinfo.randVars.percvalence(itrial); % 1 = fear, 2 = neutral
  hitmiss    = fearneut == respval; % 1 = hit, 0 = miss
  rt         = mgltrialinfo.reactionTime(itrial);
  faceindex  = mgltrialinfo.randVars.faceindex(itrial); 
  
  % get various latencies
  faceonset  = (endsample-(begsample-offset)+1) ./ hdr.Fs; % t=0 is sample=1(-offset) 
  respcueons = mgltrialinfo.trials(itrial).segtime(10)-mgltrialinfo.trials(itrial).segtime(2); % cue = t=0
  % respcueons could be on a timepoint that doesn't match a sample, find closest one
  respcueons = round(respcueons .* hdr.Fs) ./ hdr.Fs;
  
  % grab duration of cue and face+mask
  cuediodur  = event(curreventind).duration;
  facediodur = event(curreventind+1).duration;
  
  % put all in trl
  trl(end+1,:) = [begsample endsample offset predunpred fearneut hitmiss rt faceonset respcueons faceindex cuediodur facediodur];
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SUBFUNCTION
%%%%%%% this function is taken from the MGL toolbox developed by the Gardner lab at Stanford
% getTaskParameters.m
%
%      usage: exp = getTaskParameters(myscreen,task)
%         by: justin gardner
%       date: 01/27/07
%    purpose: get the task parameters, reaction times etc,
%             out of the screen and task variables
%
%             Note that all volume numbers represent the beginning of a trial or a segment
%             and are rounded to the **closest** volume number. Thus if your trial or segment
%             started at time 0.76 seconds and your frame period (TR) was 1.5 seconds, then
%             you would see a volume number of 2 rather than 1.
%
function [experiment stimfile] = getTaskParameters(myscreen,task)

if isfield(myscreen,'uniqueId')
  % this is a JGL file
  [experiment stimfile] = jglGetTaskParameters(myscreen,task);
  return
end

% check arguments
experiment = [];stimfile=[];
if ~any(nargin == [1 2])
  help getTaskParameters
  return
end

% see if we are passed the name of a file
if (nargin == 1) && isstr(myscreen)
  % check for file
  [pathStr filename ext] = fileparts(myscreen);
  if ~isempty(pathStr)
    filename = sprintf('%s.mat',fullfile(pathStr,filename));
  else
    filename = sprintf('%s.mat',fullfile(pwd,filename));
  end
  
  if ~isfile(filename)
    disp(sprintf('(getTaskParameters) Could not find file %s',filename));
    return
  end
  % load file
  load(filename);
  stimfileName = filename;
else
  stimfileName = '';
end

% so you can pass in stimfile strucuture from MLR
if isfield(myscreen,'myscreen') && isfield(myscreen,'task')
  task = myscreen.task;
  myscreen = myscreen.myscreen;
end

% get task from myscreen
if ~exist('task','var') && isfield(myscreen,'task')
  task = myscreen.task;
end

if ~exist('task','var') || isempty(task)
  disp(sprintf('(getTaskParameters) No task variable'));
  help getTaskParameters;
  return
end

% if there is only one task
if ~iscell(task)
  allTasks{1}{1} = task;
  multiTask = 0;
elseif ~iscell(task{1}) % assumes multiple phases, not tasks
  allTasks{1} = task;
  multiTask = 0;
else
  % otherwise cycle through tasks
  allTasks = task;
  multiTask = 1;
end

volumeTR = [];

for taskNum = 1:length(allTasks)
  task = allTasks{taskNum};
  % init some variables
  exptStartTime = inf;
  volnum = 0;
  nextVolNum = 1;
  volTime = 0;
  nextVolTime = inf;
  phaseNum = 1;
  blockNum = 1;
  blockTrialNum = 0;
  numTraces = max(0,max(myscreen.events.tracenum) - myscreen.stimtrace + 1);
  experiment = initPhase([],phaseNum,numTraces,task{phaseNum});
  tnum = 0;
  
  if (task{phaseNum}.segmentTrace)
    % go through the events, looking for the segment
    for enum = 1:myscreen.events.n
      % get the volume number of the event
      volnum = myscreen.events.volnum(enum);
      eventTime = myscreen.events.time(enum);
      % if we are closer to the next volume than the current
      % volume (i.e. the one last recorded by a backtick, then
      % we want to use the following volume as our volume number
      if (eventTime-volTime) > (nextVolTime-eventTime)
        volnum = nextVolNum;
      end
      % deal with segment trace
      if myscreen.events.tracenum(enum) == task{phaseNum}.segmentTrace
        if (round(myscreen.events.data(enum)) ~= myscreen.events.data(enum)) || (myscreen.events.data(enum)<1)
          disp(sprintf('(getTaskParameters) Bad segmentTrace (%i) value %i at %i',myscreen.events.tracenum(enum),myscreen.events.data(enum),enum));
          continue;
        end
        % get the segment and the segment time
        thisseg = myscreen.events.data(enum);
        segtime = myscreen.events.time(enum);
        ticknum = myscreen.events.ticknum(enum);
        % check for new trial
        if thisseg == 1
          tnum = tnum+1;
          if tnum > task{phaseNum}.numTrials
            fprintf('Recorded trace events past end of last trial.\n');
            return
          end
          experiment(phaseNum).nTrials = tnum;
          % get the time that the experiment starts
          % this will only get set for the 1st seg of 1st trial
          exptStartTime = min(segtime,exptStartTime);
          % now keep the trial time
          experiment(phaseNum).trialTime(tnum) = segtime-exptStartTime;
          experiment(phaseNum).trialTicknum(tnum) = myscreen.events.ticknum(enum);
          experiment(phaseNum).trialVolume(tnum) = volnum;
          % get block trial numbers
          blockTrialNum = blockTrialNum+1;
          % see if we have to go over to the next block
          if task{phaseNum}.block(blockNum).trialn < blockTrialNum
            blockNum = blockNum + 1;
            blockTrialNum = 1;
          end
          % save the block num and trial num
          experiment(phaseNum).blockNum(tnum) ...
            = blockNum;
          experiment(phaseNum).blockTrialNum(tnum) ...
            = blockTrialNum;
          % and initalize other parameters
          experiment(phaseNum).trials(tnum).response = [];
          experiment(phaseNum).trials(tnum).responseVolume = [];
          experiment(phaseNum).trials(tnum).responseSegnum = [];
          experiment(phaseNum).trials(tnum).reactionTime = [];
          experiment(phaseNum).trials(tnum).responseTimeRaw = [];
          experiment(phaseNum).trials(tnum).traces.tracenum = [];
          experiment(phaseNum).trials(tnum).traces.val = [];
          experiment(phaseNum).trials(tnum).traces.time = [];
          if numTraces > 0
            experiment(phaseNum).traces(:,tnum) = nan;
          end
          experiment(phaseNum).response(tnum) = nan;
          experiment(phaseNum).responseVolume(tnum) = nan;
          experiment(phaseNum).reactionTime(tnum) = nan;
          experiment(phaseNum).responseTimeRaw(tnum) = nan;
          % get all the random parameter
          for rnum = 1:task{phaseNum}.randVars.n_
            eval(sprintf('experiment(phaseNum).randVars.%s(tnum) = task{phaseNum}.randVars.%s(mod(tnum-1,task{phaseNum}.randVars.varlen_(%i))+1);',task{phaseNum}.randVars.names_{rnum},task{phaseNum}.randVars.names_{rnum},rnum));
          end
          if isfield(task{phaseNum},'parameterCode')
            experiment(phaseNum).parameterCode = task{phaseNum}.parameterCode;
          end
          % and get all parameters
          parameterNames = fieldnames(task{phaseNum}.block(blockNum).parameter);
          % and set the values
          for pnum = 1:length(parameterNames)
            thisParam = task{phaseNum}.block(blockNum).parameter.(parameterNames{pnum});
            % if it is an array then it is just a regular parameter
            if size(thisParam,1) == 1
              eval(sprintf('experiment(phaseNum).parameter.%s(tnum) = thisParam(blockTrialNum);',parameterNames{pnum}));
              % otherwise there are multiple values per each trial
            else
              for paramRowNum = 1:size(thisParam,1)
                eval(sprintf('experiment(phaseNum).parameter.%s%i(tnum) = thisParam(paramRowNum,blockTrialNum);',parameterNames{pnum},paramRowNum));
              end
            end
          end
        end
        
        % set the segment time for this trial
        segtime = segtime-exptStartTime;
        experiment(phaseNum).trials(tnum).segtime(thisseg) = segtime;
        experiment(phaseNum).trials(tnum).volnum(thisseg) = volnum;
        experiment(phaseNum).trials(tnum).ticknum(thisseg) = ticknum;
        % deal with volnum event
      elseif myscreen.events.tracenum(enum) == 1
        % if data is set to one then it means that we got a backtick
        % if it is set to zero it means we are coming out of a backtick
        if myscreen.events.data(enum)
          % remember the time of the volume
          volTime = myscreen.events.time(enum);
          % get the next volume time, by looking for the next volume event
          volEvents = find((myscreen.events.tracenum(enum+1:end) == 1) & (myscreen.events.data(enum+1:end) == 1));
          % if we have the next volume event get the time
          if ~isempty(volEvents)
            nextVolEvent = volEvents(1)+enum;
            nextVolTime = myscreen.events.time(nextVolEvent);
            nextVolNum = myscreen.events.volnum(nextVolEvent)+1;
          else
            % if we have collected some information about volumeTR
            % then we set the final+1 volume to happen one volume
            % later. This way events that happen after the last volume
            % can be set to have a volume number of nan
            if ~isempty(volumeTR(~isnan(volumeTR)))
              nextVolTime = volTime+median(volumeTR(~isnan(volumeTR)));
            else
              nextVolTime = inf;
            end
            nextVolNum = nan;
          end
          % keep the amount of time each volume takes
          volumeTR(end+1) = nextVolTime-volTime;
        end
        % deal with phasenum event
      elseif myscreen.events.tracenum(enum) == task{phaseNum}.phaseTrace
        phaseNum = myscreen.events.data(enum);
        if phaseNum <= length(task)
          blockNum = 1;
          blockTrialNum = 0;
          experiment = initPhase(experiment,phaseNum,numTraces,task{phaseNum});
          experiment(phaseNum).nTrials = 1;
          tnum = 0;
        else
          break;
        end
        % deal with response
      elseif myscreen.events.tracenum(enum) == task{phaseNum}.responseTrace
        whichButton = myscreen.events.data(enum);
        % make sure this is happening after first trial
        if tnum
          % reaction time relative to beginning of segment
          reactionTime = myscreen.events.time(enum)-exptStartTime-segtime;
          % responseTimeRaw is response time relative to beginning of the experiment
          responseTimeRaw = myscreen.events.time(enum)-exptStartTime;
          % now adjust reactionTime if the previous segments
          % had the response on (that is, the response time
          % is the time not necesarily from the beginning of
          % this current segment, but from the first segment
          % before this one in which the subject could have
          % responded.
          if isfield(task{phaseNum},'getResponse') && (length(task{phaseNum}.getResponse) >= thisseg)
            % get how long each segment took relative to the previous one
            seglen = diff(experiment(phaseNum).trials(tnum).segtime);
            % now cycle backwards from the segment previous to this one
            i = thisseg-1;
            % if getResponse was on (and we are not at the beginning of the trial yet
            while ((i >= 1) && task{phaseNum}.getResponse(i))
              % then the reaction time should be increased by the length
              % of time the segment took.
              reactionTime = reactionTime+seglen(i);
              i=i-1;
            end
          end
          % save the first response in the response array
          if isnan(experiment(phaseNum).response(tnum))
            experiment(phaseNum).response(tnum) = whichButton;
            experiment(phaseNum).reactionTime(tnum) = reactionTime;
            experiment(phaseNum).responseTimeRaw(tnum) = responseTimeRaw;
            % now see if the response happened closer to this volume
            % or closer to the next volume
            responseTime = myscreen.events.time(enum);
            experiment(phaseNum).responseVolume(tnum) = volnum;
          end
          % save all responses in trial
          experiment(phaseNum).trials(tnum).response(end+1) = whichButton;
          experiment(phaseNum).trials(tnum).reactionTime(end+1) = reactionTime;
          experiment(phaseNum).trials(tnum).responseTimeRaw(end+1) = responseTimeRaw;
          experiment(phaseNum).trials(tnum).responseSegnum(end+1) = thisseg;
          experiment(phaseNum).trials(tnum).responseVolume(end+1) = volnum;
        end
        % deal with user traces
      elseif myscreen.events.tracenum(enum) >= myscreen.stimtrace
        tracenum = myscreen.events.tracenum(enum)-myscreen.stimtrace+1;
        userval = myscreen.events.data(enum);
        usertime = myscreen.events.time(enum)-exptStartTime;
        % there is some chance that a user trace can be written
        % before the first trial is started for this task. This
        % happens if there are multiple tasks and this user
        % trace belongs to another task. In that case, storing
        % this variable with this task is not really necessary,
        % but we do not know that here so we just either save
        % it if we have a valid trial number or ignore it if not.
        if (tnum)
          % store it if it is the first setting
          if isnan(experiment(phaseNum).traces(tracenum,tnum))
            experiment(phaseNum).traces(tracenum,tnum) = userval;
          end
          % put it in trial
          experiment(phaseNum).trials(tnum).traces.tracenum(end+1) = tracenum;
          experiment(phaseNum).trials(tnum).traces.val(end+1) = userval;
          experiment(phaseNum).trials(tnum).traces.time(end+1) = usertime;
        end
      end
    end
  end
  % for a multi task experiment, then we keep a cell array of values
  if multiTask
    retval{taskNum} = experiment;
  else
    retval = experiment;
  end
end


experiment = retval;

% save stimfile
stimfile.stimfilePath = '';
if ~isempty(stimfileName)
  [stimfile.stimfilePath stimfile.stimfile] = fileparts(stimfileName);
else
  if isfield(myscreen,'stimfile');
    stimfile.stimfile = myscreen.stimfile;
  else
    stimfile.stimfile = '';
  end
end
stimfile.myscreen = myscreen;
stimfile.task = allTasks;

% set the traces in the return value if they exist
if isfield(myscreen,'traces')
  if iscell(experiment)
    for i = 1:length(experiment)
      for j = 1:length(experiment{i})
        experiment{i}(j).tracesAll = myscreen.traces;
      end
    end
  else
    for j = 1:length(experiment)
      experiment(j).tracesAll = myscreen.traces;
    end
  end
end

%%%%%%%%%%%%%%%%%%%
%    initPhase    %
%%%%%%%%%%%%%%%%%%%
function experiment = initPhase(experiment,phaseNum,numTraces,task)

experiment(phaseNum).nTrials = 0;
experiment(phaseNum).trialVolume = [];
experiment(phaseNum).trialTime = [];
experiment(phaseNum).trialTicknum = [];
experiment(phaseNum).trials = [];
experiment(phaseNum).blockNum = [];
experiment(phaseNum).blockTrialNum = [];
experiment(phaseNum).response = [];
experiment(phaseNum).reactionTime = [];
if numTraces>0
  experiment(phaseNum).traces(1:numTraces,:) = nan;
end

% get what the parameters were originaly set to - i.e. in the task variable. This gives a record
% of all the values that the parameter was originally intended to go through for example (sometimes
% an experiment may run short and you don't go through all possible values).
experiment(phaseNum).originalTaskParameter = task.parameter;
taskParameters = fieldnames(task.parameter);
for i = 1:length(taskParameters)
  % remove fields that end in _ which are created by initRandomization
  if taskParameters{i}(end) == '_'
    experiment(phaseNum).originalTaskParameter = rmfield(experiment(phaseNum).originalTaskParameter,taskParameters{i});
  else
    % check for a multi-row field, this a variable that has different settings for each row
    % i.e. like when you do a split screen design with one randomization for the left and one for the right
    numRows = size(task.parameter.(taskParameters{i}),1);
    if numRows > 1
      % remove the field
      experiment(phaseNum).originalTaskParameter = rmfield(experiment(phaseNum).originalTaskParameter,taskParameters{i});
      % and reset to a field name which has a number for each row. e.g. orientation will become
      % orientation1, orientation2 etc.
      for iRows = 1:numRows
        experiment(phaseNum).originalTaskParameter.(sprintf('%s%i',taskParameters{i},iRows)) = task.parameter.(taskParameters{i})(iRows,:);
      end
    end
  end
end

% now do the same for randVars. Works the same way, but get what the parameters were originaly set to - i.e. in the task variable. This gives a record
randVarTypes = {'uniform','block','calculated'};
for iRandVarType = 1:length(randVarTypes)
  if isfield(task.randVars,randVarTypes{iRandVarType})
    randVars = fieldnames(task.randVars.(randVarTypes{iRandVarType}));
    for i = 1:length(randVars)
      % only use fields that don't end in _
      if randVars{i}(end) ~= '_'
        % check to see if there is a variable with the same name except with an underscore after
        % it, that will contain the variables all possible settings.
        allSettings = find(strcmp(sprintf('%s_',randVars{i}),randVars));
        if ~isempty(allSettings)
          experiment(phaseNum).originalRandVars.(randVars{i}) = task.randVars.(randVarTypes{iRandVarType}).(randVars{allSettings});
        else
          % otherwise just set to whatever it was set to
          experiment(phaseNum).originalRandVars.(randVars{i}) = task.randVars.(randVarTypes{iRandVarType}).(randVars{i});
        end
      end
    end
  end
end
























