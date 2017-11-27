function data = rmr_crcnsssc3_readspike(sessionpath,triallength)

% This function reads and processes the spike data from CRCNS SSC3 datasets
%
% It results in a FieldTrip style data structure (Note, not a spike data structure)
% Where the spike times are represented as a boolean spike train. Each unit is
% a separate channel.
%
%
% The data is segmented into trials.
% Two exceptions w.r.t. regular fieldtrip data structure. 
% 1) data.trial is sparse
% 2) data.time only contains [begtime endtime]
%
% Spike data is sampled at 20kHz
%
% %%%
% CRUCIAL(!): time stamps of events and spikes MUST be coded relative to the start of the experiment 
% %%%
%
%        % trialinfo variables:
%            1st column: trial number
%
%
%



% set paths/filenames and load data
%sessionpath = '/Users/roemervandermeij/Work/Data/tmpdata/CRCNS_SSC3/DataSet1.mat';
rawdata = load([sessionpath '.mat']);
rawdata = rawdata.data;

% set basics
nsample    = rawdata.recordinglength*20;
fsample    = 20000;  %%% HARDCODED
nunit      = rawdata.nNeurons;
datalength = nsample ./ fsample; % in seconds


%%%%%%%%%% Create trl
% create trial boundaries
begsample = 1:round(triallength.*fsample):nsample;
endsample = (round(triallength.*fsample)):round(triallength.*fsample):nsample;
begsample = begsample(1:min(numel(begsample),numel(endsample)));
endsample = endsample(1:min(numel(begsample),numel(endsample)));

% get trl
trl = [begsample' endsample' zeros(numel(begsample),1) (1:numel(begsample))'];

ntrial = size(trl,1);
%%%%%%%%%% Create trl


%%%%%%%%%% Create spiketrains 
% create trial and time cell-array containing spike trains
trial   = cell(1,ntrial);
time    = cell(1,ntrial);
for itrial = 1:ntrial
 
  % create trial-specific time axis
  begsample = trl(itrial,1);
  endsample = trl(itrial,2);
  nsample   = (endsample-begsample)+1;
  currtime  = ([begsample endsample]-begsample) ./ fsample; % HUGE EXCEPTION
  %currtime  = (1:nsample) ./ fsample; 
     
  % for each unit, transform the spike timestamp samples into a sparsified spike train
  spunitind = []; % unit index for creating sparse trial matrix
  sptsind   = []; % time point/stamp index for creating sparse trial matrix
  spvalind  = []; % values for creating sparse trial matrix (all 1s of course)
  for iunit = 1:nunit
    
    % get timestamp samples for current unit 
    currtss = rawdata.spikes{iunit}*20;
    if ~isequal(round(currtss),double(int64(currtss)))
      error('woopsie')
    else
      currtss = round(currtss);
    end
    currtss = currtss(currtss >= begsample & currtss <= endsample);
    currtss = (currtss - begsample) + 1; % convert to trial specific samples (i.e. indices)
    if ~isempty(currtss)
      % add to sparsified spike trains
      nspikes = numel(currtss);
      spunitind = [spunitind ones(1,nspikes) .* iunit];
      sptsind   = [sptsind currtss(:)'];
      spvalind  = [spvalind ones(1,nspikes)];
    end
  end
  % create sparse representation of spiketrains
  currtrial = sparse(spunitind,sptsind,spvalind,nunit,nsample);
  % save
  trial{itrial} = currtrial;
  time{itrial}  = currtime;
end
%%%%%%%%%%


%%%%%%%%%% Create label field
label = [];
for iunit = 1:nunit
  label{iunit} = ['unit' num2str(iunit) '_' 'x' num2str(round(rawdata.x(iunit))) 'y' num2str(round(rawdata.y(iunit)))];
end
%%%%%%%%%%



%%%%%%%%%% Create data structure
data            = [];
data.hdr.Fs     = fsample;
data.label      = label;
data.time       = time;
data.trial      = trial;
data.fsample    = fsample;
data.sampleinfo = round(trl(:,1:2));
data.trialinfo  = trl(:,4:end);
data.cfg        = [];
data.cfg.trl    = trl;
%%%%%%%%%%




























