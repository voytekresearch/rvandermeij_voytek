function data = rmr_carmonkey_readspike(sessionpath,fsample,trialdeftype)

% This function reads and processes the spike data from monkey Paco from the Carmena lab
%
% It results in a FieldTrip style data structure (Note, not a spike data structure)
% Where the spike times are represented as a boolean spike train. Each unit is
% a separate channel.
%
% The data is segmented into trials.
% Two exceptions w.r.t. regular fieldtrip data structure.
% 1) data.trial is sparse
% 2) data.time only contains [begtime endtime]
%
% Spike data was originally sampled at 40kHz
%
% %%%
% CRUCIAL(!): time stamps of events and spikes MUST be coded relative to the start of the experiment
% %%%
%
% From the loaded the data, the following variables contain:
%   AD39-AD40 = cursor coordinates (X,Y?) (not necessary now)
%    AD65-192 = LFP channels 1:128
%     Strobed =  2xNevent vector of events and their time-points
%
%
% W.r.t. to LFP channel identity, the following was provided by the Carmena lab:
% Channels 1-64 correspond to ipsi (R) hemisphere (except, we suspect that channels 45,15,53,55 and 59 are LPMd), and
% 65-128 are contra (L) hemisphere.
%
%
%

% error checking
if fsample<40000
  error('sampling rate lower than spike sampling rate not supported')
end

% set paths and filenames
% datapath = '/Users/roemer/Work/Data/carmena_monkey/';
% sessionname = 'paco020608c'; % 'paco020608b'

% load data
load([sessionpath '_AD33_AD40.mat'])
load([sessionpath '_AD65_AD128.mat'])
load([sessionpath '_AD129_AD192.mat'])
load([sessionpath '_events.mat'])
load([sessionpath '_sig.mat'])
load([sessionpath '_signal_parameters.mat'])



%%%%%%%%%% Create trl
% define trials using the LFP sampling
% define trial starting with event START (2), and ends with event SUCCESS (11), TIME_OUT (12) or TARGET_HOLD_ERROR (8)
% after a START event a BC_TX (73-80) event follows, indicating which location is the target
trl = [];
eventcode = Strobed(:,2);
eventtime = Strobed(:,1);
nevent    = numel(eventcode);
% loop over events, and create trl matrix
count = 0;
for ievent = 1:nevent
  
  switch trialdeftype
    
    case 'tarbackw'
      % hook to START(2) event
      if eventcode(ievent)==2
        startievent = ievent;
        % search for end event, being either SUCCESS (11), TIME_OUT (12) or TARGET_HOLD_ERROR (8)
        endfound  = false;
        currievent = ievent + 1;
        while ~endfound && currievent<=nevent && eventcode(currievent)~=2 % stop search if another START event is detected
          if any(ismember(eventcode(currievent),[11])) % 11 is success, 12 is timeout, 8 is targethold error
            % succes!
            endievent = currievent;
            endfound  = true;
            
            % set endievent to either 7 (entering target) or 12
            if eventcode(endievent)==11
              endievent = startievent+5; % should always be 7
              if eventcode(endievent)~=7
                error('event is not 7!')
              end
            elseif eventcode(endievent)==12
              % do nothing
            elseif eventcode(endievent)==8
              endievent = startievent+5; % should always be 7
              if eventcode(endievent)~=7
                error('event is not 7!')
              end
            end
            
            % define the trial - going backwards from entering target
            % the start sample, end sample and offset
            begsample = round(eventtime(endievent-3) .* fsample)+1; % EVENT 5, (in case of SUCCESS)
            endsample = round(eventtime(endievent-1) .* fsample)+1; % EVENT 7, (in case of SUCCESS)
            begsample = max([begsample endsample-round(0.5*fsample)]); % cap at 500ms seconds
            offset = 0;
            
            % event codes for the key events
            finalevent = eventcode(endievent);
            startevent = eventcode(endievent-3);
            endevent   = eventcode(endievent-1);
            
            % trailseq, the sequence of events of the trial
            trialseq = eventcode(startievent:endievent)';
            
            % get TIMETILLNEXTEVENT in ms
            trialtim = round(diff(eventtime(startievent:(endievent+1))).*1000);
            
            % merge trialtim with trialseq to display next to each other
            trialst = [trialseq' trialtim]';
            trialst = trialst(:)';
            
            % put in trl
            trl(end+1,1:(numel(trialseq)*2)+6) = [begsample endsample offset finalevent startevent endevent trialst];
          else
            currievent = currievent + 1;
          end
        end
      end
      
      
    case 'gotilltar'
      % hook to START(2) event
      if eventcode(ievent)==2
        startievent = ievent;
        % search for end event, being either SUCCESS (11), TIME_OUT (12) or TARGET_HOLD_ERROR (8)
        endfound  = false;
        currievent = ievent + 1;
        while ~endfound && currievent<=nevent && eventcode(currievent)~=2 % stop search if another START event is detected
          if any(ismember(eventcode(currievent),[11 12 8])) % 11 is success, 12 is timeout, 8 is targethold error
            % succes!
            endievent = currievent;
            endfound  = true;
            
            % set endievent to either 7 (entering target) or 12 (timeout)
            if eventcode(endievent)==11
              endievent = startievent+5; % should always be 7
              if eventcode(endievent)~=7
                error('event is not 7!')
              end
            elseif eventcode(endievent)==12
              % do nothing
            elseif eventcode(endievent)==8
              endievent = startievent+5; % should always be 7
              if eventcode(endievent)~=7
                error('event is not 7!')
              end
            end
            
            % define the trial
            % the start sample, end sample and offset
            begsample = round(eventtime(startievent+2) .* fsample)+1; % EVENT 15, always comes two events later than the start
            endsample = round(eventtime(endievent)     .* fsample)+1;
            endsample = min([endsample round(eventtime(endievent-1) .* fsample)+1+(10*fsample)]); % cap at 10 seconds
            offset = 0;
            
            % event codes for the key events
            finalevent = eventcode(endievent);
            startevent = eventcode(startievent+2);
            endevent   = eventcode(endievent);
            
            % trailseq, the sequence of events of the trial
            trialseq = eventcode(startievent:endievent)';
            
            % get TIMETILLNEXTEVENT in ms
            trialtim = round(diff(eventtime(startievent:(endievent+1))).*1000);
            
            % merge trialtim with trialseq to display next to each other
            trialst = [trialseq' trialtim]';
            trialst = trialst(:)';
            
            % put in trl
            trl(end+1,1:(numel(trialseq)*2)+6) = [begsample endsample offset finalevent startevent endevent trialst];
          else
            currievent = currievent + 1;
          end
        end
      end
      
      
    case 'gotillend'
      % hook to START(2) event
      if eventcode(ievent)==2
        startievent = ievent;
        % search for end event, being either SUCCESS (11), TIME_OUT (12) or TARGET_HOLD_ERROR (8)
        endfound  = false;
        currievent = ievent + 1;
        while ~endfound && currievent<=nevent && eventcode(currievent)~=2 % stop search if another START event is detected
          if any(ismember(eventcode(currievent),[11 12 8])) % 11 is success, 12 is timeout, 8 is targethold error
            % succes!
            endievent = currievent;
            endfound  = true;
            
            % define the trial
            % the start sample, end sample and offset
            begsample = round(eventtime(startievent+2) .* fsample)+1; % EVENT 15, always comes two events later than the start
            endsample = round(eventtime(endievent)     .* fsample)+1;
            endsample = min([endsample round(eventtime(endievent-1) .* fsample)+1+(10*fsample)]); % cap at 10 seconds
            offset = 0;
            
            % event codes for the key events
            finalevent = eventcode(endievent);
            startevent = eventcode(startievent+2);
            endevent   = eventcode(endievent);
            
            % trailseq, the sequence of events of the trial
            trialseq = eventcode(startievent:endievent)';
            
            % get TIMETILLNEXTEVENT in ms
            trialtim = round(diff(eventtime(startievent:(endievent+1))).*1000);
            
            % merge trialtim with trialseq to display next to each other
            trialst = [trialseq' trialtim]';
            trialst = trialst(:)';
            
            % put in trl
            trl(end+1,1:(numel(trialseq)*2)+6) = [begsample endsample offset finalevent startevent endevent trialst];
          else
            currievent = currievent + 1;
          end
        end
      end
      
      
    case 'starttilltar'
      
      % hook to START(2) event
      if eventcode(ievent)==2
        startievent = ievent;
        % search for end event, being either SUCCESS (11), TIME_OUT (12) or TARGET_HOLD_ERROR (8)
        endfound  = false;
        currievent = ievent + 1;
        while ~endfound && currievent<=nevent && eventcode(currievent)~=2 % stop search if another START event is detected
          if any(ismember(eventcode(currievent),[11 12 8])) % 11 is success, 12 is timeout, 8 is targethold error
            % succes!
            endievent = currievent;
            endfound  = true;
            
            % set endievent to either 7 (entering target) or 12
            if eventcode(endievent)==11
              endievent = startievent+5; % should always be 7
              if eventcode(endievent)~=7
                error('event is not 7!')
              end
            elseif eventcode(endievent)==12
              % do nothing
            elseif eventcode(endievent)==8
              endievent = startievent+5; % should always be 7
              if eventcode(endievent)~=7
                error('event is not 7!')
              end
            end
            
            % define the trial
            % the start sample, end sample and offset
            begsample = round(eventtime(startievent)  .* fsample)+1;
            endsample = round(eventtime(endievent)    .* fsample)+1;
            offset = 0;
            
            % event codes for the key events
            finalevent = eventcode(endievent);
            startevent = eventcode(startievent);
            endevent   = eventcode(endievent);
            
            % trailseq, the sequence of events of the trial
            trialseq = eventcode(startievent:endievent)';
            
            % get TIMETILLNEXTEVENT in ms
            trialtim = round(diff(eventtime(startievent:(endievent+1))).*1000);
            
            % merge trialtim with trialseq to display next to each other
            trialst = [trialseq' trialtim]';
            trialst = trialst(:)';
            
            % put in trl
            trl(end+1,1:(numel(trialseq)*2)+6) = [begsample endsample offset finalevent startevent endevent trialst];
          else
            currievent = currievent + 1;
          end
        end
      end
      
      
    case 'centertillgo'
      
      % hook to START(2) event
      if eventcode(ievent)==2
        startievent = ievent;
        % search for end event, being either SUCCESS (11), TIME_OUT (12) or TARGET_HOLD_ERROR (8)
        endfound  = false;
        currievent = ievent + 1;
        while ~endfound && currievent<=nevent && eventcode(currievent)~=2 % stop search if another START event is detected
          if any(ismember(eventcode(currievent),[11 12 8])) % 11 is success, 12 is timeout, 8 is targethold error
            % succes!
            endievent = currievent;
            endfound  = true;
           
            % define the trial
            % the start sample, end sample and offset
            begsample = round(eventtime(startievent+2)  .* fsample)+1; % startievent+2 is always 15 (enters center)
            endsample = round(eventtime(startievent+3)  .* fsample)+1; % startievent+2 is always 5 (go cue)
            offset = 0;
            
            % event codes for the key events
            finalevent = eventcode(endievent);
            startevent = eventcode(startievent+2);
            endevent   = eventcode(startievent+3);
            
            % trailseq, the sequence of events of the trial
            trialseq = eventcode(startievent:endievent)';
            
            % get TIMETILLNEXTEVENT in ms
            trialtim = round(diff(eventtime(startievent:(endievent+1))).*1000);
            
            % merge trialtim with trialseq to display next to each other
            trialst = [trialseq' trialtim]';
            trialst = trialst(:)';
            
            % put in trl
            trl(end+1,1:(numel(trialseq)*2)+6) = [begsample endsample offset finalevent startevent endevent trialst];
          else
            currievent = currievent + 1;
          end
        end
      end
      
      
    case 'starttillendseg'
      
      % hook to START(2) event
      if eventcode(ievent)==2
        startievent = ievent;
        % search for end event, being either SUCCESS (11), TIME_OUT (12) or TARGET_HOLD_ERROR (8)
        endfound  = false;
        currievent = ievent + 1;
        while ~endfound && currievent<=nevent && eventcode(currievent)~=2 % stop search if another START event is detected
          if any(ismember(eventcode(currievent),[11 12 8])) % 11 is success, 12 is timeout, 8 is targethold error
            % succes!
            endievent = currievent;
            endfound  = true;
            count     = count + 1;
            
            % divide the trial into individual segments
            segevents = startievent:endievent;
            segevents(2) = []; % always ditch the first event after start, target code, which monkey doesn't see
            segevents(eventcode(segevents)==51) = []; % ditch the mystery events
            segbegievent = segevents(1:end-1)';
            segendievent = segevents(2:end)';
            % throw away really short segments, i.e. 10ms
            remind = (eventtime(segendievent) - eventtime(segbegievent)) <0.010;
            segbegievent(remind) = [];
            segendievent(remind) = [];
            nseg = numel(segbegievent);
            
            % first 3 columns of trl
            begsample = round(eventtime(segbegievent)  .* fsample)+1;
            endsample = round(eventtime(segendievent)  .* fsample);
            offset    = ones(nseg,1) .* 0;
            % cap start-to-enters-center segment (always 1st) at 15s, same as timeout period, going backwards
            begsample(1) = max(begsample(1),endsample(1)-(15*fsample));
            
            % trial id
            trialnum = ones(nseg,1) .* count;
            
            % event code for the last event
            endevent = ones(nseg,1) .* eventcode(endievent);
            
            % start/end event codes for the segment
            segbegevent = eventcode(segbegievent);
            segendevent = eventcode(segendievent);
            
            % target location
            targetloc = ones(nseg,1) .* eventcode(startievent+1);
            
            % put in trl
            trl(end+1:end+nseg,1:8) = [begsample endsample offset trialnum endevent segbegevent segendevent targetloc];
          else
            currievent = currievent + 1;
          end
        end
      end
      
      
    otherwise
      error('trialdeftype not supported')
  end
end
%%%%%%%%%% Create trl



%%%%%%%%%% Create spiketrains
% create trial and time cell-array containing spike trains
ntrial = size(trl,1);
filevars  = whos('-file', [sessionpath '_sig.mat']);
label = {filevars.name};
label = label(strncmp('sig',label,3));
% remove sigXXXi (unsorted units)
remind = [];
for ilab = 1:numel(label)
  if strcmp(label{ilab}(end),'i')
    remind = [remind ilab];
  end
end
label(remind) = [];
nunit = numel(label);
%
trial = cell(1,ntrial);
time  = cell(1,ntrial);
for itrial = 1:ntrial
  % create trial-specific time axis
  begsample = trl(itrial,1);
  endsample = trl(itrial,2);
  currtime  = ([begsample endsample]-begsample) ./ fsample; % HUGE EXCEPTION
  nsample   = (endsample-begsample)+1;
  % for each unit, transform the spike timestamps into a sparsified spike train
  spunitind = []; % unit index for creating sparse trial matrix
  sptsind   = []; % time point/stamp index for creating sparse trial matrix
  spvalind  = []; % values for creating sparse trial matrix (all 1s of course)
  for iunit = 1:nunit
    eval(['currts = ' label{iunit} ';']); % that it has come to this...
    currts = round(currts .* fsample) + 1; % convert to samples
    ind    = currts>=begsample & currts<=endsample;
    currts = (currts(ind) - begsample) + 1; % convert to trial specific samples (i.e. indices)
    if ~isempty(currts)
      % add to sparsified spike trains
      nspikes = numel(currts);
      spunitind = [spunitind ones(1,nspikes) .* iunit];
      sptsind   = [sptsind currts(:)'];
      spvalind  = [spvalind ones(1,nspikes)];
    end
  end
  % create sparse representation of spiketrains
  currtrial = sparse(spunitind,sptsind,spvalind,nunit,nsample);
  % save
  trial{itrial} = currtrial;
  time{itrial}  = currtime;
end
%%%%%%%%%% Create spiketrains



%%%%%%%%%% Create data structure
data = [];
data.hdr.Fs     = fsample;
data.label      = label;
data.time       = time;
data.trial      = trial;
data.fsample    = fsample;
data.sampleinfo = trl(:,1:2);
data.trialinfo  = trl(:,4:end);
data.cfg        = [];
data.cfg.trl    = trl;
data.cfg.previous.trialdeftype = trialdeftype;
data.cfg.previous.fsample      = fsample;
%%%%%%%%%%









