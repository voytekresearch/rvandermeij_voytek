
function data = rmr_crcnshc5_readspike(sessionpath,fsample,trialdeftype,plotresult,incintern,nsegments)

% This function reads and processes the spike data from CRCNS HC5 datasets
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
% Currently supported trial definitions (given as a string in trialdeftype:
%   - *mazerun* 
%     Each lap (nlap = ~15) is segemented into several sections based on the linearized DISTANCE of the 
%     rat along the track (not the time during a lap). A typical lap takes about 4-5 sec. Distance along the track is divided
%     into segments 10% (of total distance) with 50% overlap
%     and then segmented with 50% overlap in .
%        % trialinfo variables:
%            1st column: lap number
%            2st column: lap segment number
%            3rd column: 1 (left) or 2 (right) lap (after the ~8th segment the rat is in either the L or R part)
%            4th column: Par.BehavType for this lap (code unknown)
%            5th column: Par.TrialType for this lap (code unknown)
%   
%   
%   
%   

% input and error checking
if fsample<20000
  error('sampling rate lower than spike sampling rate not supported')
end
plotresult = istrue(plotresult);
incintern  = istrue(incintern);

% set paths/filenames and load data
%sessionpath = '/Users/roemer/Work/Data/CRCNS_HC5/i01_maze06.002/';
% find last part of fn
if ~strcmp(sessionpath(end),'/')
  sessionpath = [sessionpath '/'];
end
indsl = strfind(sessionpath,'/');
inddt = strfind(sessionpath,'.');
datfn = sessionpath(indsl(end-1)+1:inddt(end)-1);
datfn = [datfn '_MS' sessionpath(inddt(end):end-1) '_BehavElectrData.mat']; 
load([sessionpath datfn])


% get fsamplefac for transforming 'undersampled' events and spike time stamps to target sampling rate
eventfsfac = fsample ./ Par.SamplingFrequency;
spikefsfac = fsample ./ 20000; % sampling rate hardcoded into Spike.res20khz



%%%%%%%%%% Create trl 
% define trials using Par, Track, Laps
% specify trl in time stamps at desired resolution
trl = [];
switch trialdeftype
 
  case 'mazerun'
    nlap    = Par.NLap;
    nsample = Par.nTimebins;
    
    %%%%%%%% MAX LINEARIZED DISTANCE IS HARDCODED
    maxdist = 2300; % threshold considered 'end of track' (hardcoded, cause currently unknown where lindist comes from)
    mindist = 100;
    %%%%%%%%
    
    % loop over laps, further segment using Track.lindist and create trl matrix
    for ilap = 1:nlap
      % get lap start (this usually around lindist = 100)
      beglapsmp = Par.MazeSectEnterLeft{ilap}(1,1);
      % get lap end (search capped at 50 sec, 10x more than expected)
      distfromstart = Track.linDist(beglapsmp:min(beglapsmp+round(Par.SamplingFrequency*50),nsample));
      endlapsmp = beglapsmp -1 + find(distfromstart>maxdist,1); % first time point where rat is considered at 'end of track'
      
      
      % create segments
      segdist   = mindist:(maxdist-mindist)/(nsegments+1):maxdist;
      begsegdist = segdist(1:1:end-2)';
      endsegdist = segdist(3:1:end)';
      lapdist = Track.linDist(beglapsmp:endlapsmp);
      % remove segments not present
      remind = endsegdist<=lapdist(1) | begsegdist>=lapdist(end);
      begsegdist(remind) = [];
      endsegdist(remind) = [];
      begsegdist(begsegdist<lapdist(1))   = lapdist(1); % correct for beginning and end differences in lapdist
      endsegdist(endsegdist>lapdist(end)) = lapdist(end);
      nseg = numel(endsegdist);
      begsegsmp = NaN(nseg,1);
      endsegsmp = NaN(nseg,1);
      for iseg = 1:nseg
        begsegsmp(iseg) = beglapsmp -1 + find(lapdist>begsegdist(iseg),1);
        endsegsmp(iseg) = beglapsmp -1 + find(lapdist>=endsegdist(iseg),1);
      end
      
      % determine L or R movement
      avgypos = mean(Track.Y(begsegsmp(1):endsegsmp(end)));
      if avgypos<500
        LorR = 1; % rat was taking a left lap
      else
        LorR = 2; % rat was taking a right lap
      end
      
      % create laptrl
      begsample = begsegsmp;
      endsample = endsegsmp;
      offset    = zeros(nseg,1);
      lapnum    = ones(nseg,1) .* ilap;
      segnum    = (1:nseg).';
      segtime   = ((endsegsmp-begsegsmp) .* eventfsfac) ./ fsample;
      lrlap     = ones(nseg,1) .* LorR;
      behavtype = ones(nseg,1) .* Par.BehavType(ilap);
      trialtype = ones(nseg,1) .* Par.TrialType(ilap);
      laptrl = [begsample endsample offset lapnum segnum lrlap behavtype trialtype segtime begsegdist endsegdist];
   
      % add to trl
      trl(end+1:end+nseg,1:size(laptrl,2)) = laptrl;
    end
    
    
    
  otherwise
    error('undefined trialdeftype')
end % switch
ntrial = size(trl,1);
%%%%%%%%%% Create trl



%%%% Debug plotting of trial definition
if plotresult
  figure('numbertitle','off','name','trial definition result');
  plot((1:nsample)./Par.SamplingFrequency,Track.linDist);
  x = [trl(:,1) trl(:,1) trl(:,2) trl(:,2)]';
  x = x ./ Par.SamplingFrequency;
  y = trl(:,5);
  y = (y ./ max(y) .* (maxdist-100))+100;
  y = [zeros(ntrial,1) y y zeros(ntrial,1)]';
  line(x,y,'color','r');
  legend('linearized distance','trials, start to end')
  xlabel('time (s)')
  ylabel('linearized distance (au)')
end
%%%% Debug plotting of trial definition



%%%%%%%%%% Create spiketrains 
% create trial and time cell-array containing spike trains
% select only pyramidal cells
if incintern
  unitsel = 1:numel(Clu.isIntern);
  [dum sortind] = sort(Clu.isIntern); % sort according to cell type, P first, I last
  unitsel = unitsel(sortind);
else
  unitsel = find(~Clu.isIntern);
end
nunit   = numel(unitsel);
trial   = cell(1,ntrial);
time    = cell(1,ntrial);
% get timestamp samples at new sampling rate 
tsssnewfs = round(Spike.res20kHz .* spikefsfac);
for itrial = 1:ntrial
  % create trial-specific time axis
  begsample = round(trl(itrial,1) .* eventfsfac);
  endsample = round(trl(itrial,2) .* eventfsfac);
  nsample   = (endsample-begsample)+1;
  currtime  = ([begsample endsample]-begsample) ./ fsample; % HUGE EXCEPTION
  %currtime  = (1:nsample) ./ fsample; 
  
  
  % fetch timestamp samples from all units in current trial
  ind = tsssnewfs >= begsample & tsssnewfs <= endsample;
  currtrialtss = tsssnewfs(ind);
  currtrialuid = Spike.totclu(ind);
  
  % for each unit, transform the spike timestamp samples into a sparsified spike train
  spunitind = []; % unit index for creating sparse trial matrix
  sptsind   = []; % time point/stamp index for creating sparse trial matrix
  spvalind  = []; % values for creating sparse trial matrix (all 1s of course)
  for iunit = 1:nunit
    % set unit id
    curruid = unitsel(iunit);
    
    % get timestamp samples for current unit 
    currtss = currtrialtss(currtrialuid == curruid); % select unit
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
  curruid = unitsel(iunit);
  if Clu.isIntern(curruid)
    neurontype = 'I';
  else
    neurontype = 'P';
  end
  label{iunit} = ['unit' num2str(curruid) neurontype '_' num2str(Clu.shank(curruid)) '.' num2str(Clu.SpatLocalChan(curruid)) '.' num2str(Clu.localClu(curruid))];
end
%%%%%%%%%%



%%%%%%%%%% Create data structure
data            = [];
data.hdr.Fs     = fsample;
data.label      = label;
data.time       = time;
data.trial      = trial;
data.fsample    = fsample;
data.sampleinfo = round(trl(:,1:2) .* eventfsfac);
data.trialinfo  = trl(:,4:end);
data.cfg        = [];
data.cfg.trl    = trl;
%%%%%%%%%%





















function playground

% plot per trail the x-y position of the animal to determine what the hell codes for L/R laps
figure
nlap = numel(unique(trl(:,4)));
for ilap = 1:nlap
  subplot(ceil(sqrt(nlap)),ceil(sqrt(nlap)),ilap)
  % find start/end
  begsmp = trl(find(trl(:,4)==ilap,1),1);
  endsmp = trl(find(trl(:,4)==ilap,1,'last'),2);
  X = Track.X(begsmp:endsmp);
  Y = Track.Y(begsmp:endsmp);
  lindist = Track.linDist(begsmp:endsmp);
  line(Y,X)
  title(['lap ' num2str(ilap) ' LR ' num2str(trl(find(trl(:,4)==ilap,1),6))])
  axis([0 1000 0 1000])
  xlabel('Y')
  ylabel('X')
end


% plot per trail the x-y position of the animal to determine what the hell codes for L/R laps
figure
nlap = numel(unique(trl(:,4)));
for ilap = 1:nlap
  subplot(ceil(sqrt(nlap)),ceil(sqrt(nlap)),ilap)
  % find start/end
  begsmp = trl(find(trl(:,4)==ilap,1),1);
  endsmp = trl(find(trl(:,4)==ilap,1,'last'),2);
  X = Track.X(begsmp:endsmp);
  Y = Track.Y(begsmp:endsmp);
  lindist = Track.linDist(begsmp:endsmp);
  hold on
  plot(X,'b')
  plot(Y,'r')
  plot(lindist,'gr')
  title(['lap ' num2str(ilap) ' LR ' num2str(trl(find(trl(:,4)==ilap,1),6))])
  %axis([0 1000 0 1000])
end




% plot ISI's
unitsel = 1:numel(Clu.isIntern);
[dum sortind] = sort(Clu.isIntern); % sort according to cell type, P first, I last
unitsel = unitsel(sortind);
nunit   = numel(unitsel);
timbins = (0:0.0001:0.01);
count   = zeros(numel(timbins),nunit);
for iunit = 1:nunit
  curruid = unitsel(iunit);
  currspikes = sort(Spike.res20kHz(Spike.totclu==curruid) ./ 20000);
  currspikes = diff(currspikes);
  count(:,iunit) = (histc(currspikes,timbins) ./ sum(Spike.totclu==curruid))*100;
end
for iset = 1:5
  figure
  for iunit = (1:25)+((iset-1)*25)
    if iunit>nunit
      continue
    end
    subplot(5,5,iunit-((iset-1)*25))
    currcount = count(:,iunit);
    bar(timbins+diff(timbins(1:2)/2),currcount,1)
    title(label{iunit},'interpreter','none')
  end
end

figure;
bar(timbins+diff(timbins(1:2)/2),mean(count,2),1)

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



