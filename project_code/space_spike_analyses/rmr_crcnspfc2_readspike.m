function data = rmr_crcnspfc2_readspike(sessionpath,fsample,trialdeftype,plotresult,nsegments)

% This function reads and processes the spike data from CRCNS PFC2 datasets
% (only those from odor working memory maze alternation)
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
%     Each lap (nlap = 15-40) is segmented into several sections based on the linearized DISTANCE of the
%     rat along the track (not the time during a lap). A typical lap takes about 4-5 sec. Distance along the track is divided
%     into segments 10% (of total distance) with 50% overlap
%     and then segmented with 50% overlap
%        - trialinfo variables:
%            1st column: lap number
%            2st column: lap segment number
%            3rd column: 1 (right) or 2 (left) lap
%            4th column: segment length in seconds
%            5th column: position on track at start of segment
%            6th column: position on track at end of segment
%
%
%
%

% input and error checking
if fsample<20000
  error('sampling rate lower than spike sampling rate not supported')
end

% get fsamplefac for transforming 'undersampled' events and spike time stamps to target sampling rate
fsframe = 20000./512; % from docs
fsspike = 20000; % from docs
eventfsfac = fsample ./ fsframe; % in this case, events are taken from the whl files, which is per video frame (at 20k/512Hz)
spikefsfac = fsample ./ fsspike; %


% set paths/filenames and load data
% find last part of fn
if ~strcmp(sessionpath(end),'/')
  sessionpath = [sessionpath '/'];
end
ind      = strfind(sessionpath,'/');
sessname = sessionpath(ind(end-1)+1:ind(end)-1);
% treat EE.188 specially
if ~strcmp(sessname,'EE.188')
  fn = [sessionpath sessname '_Behavior.mat'];
  dat = load(fn);
else
  fn = [sessionpath 'EE188' '_example.mat'];
  dat = load(fn);
  % estimate variables to be used for trial definition
  % SessionNP
  % whlrld(:,[1 2 7])
  % Cells (~CellsInfo)
  % spikeph (~CellsInfo)
  
  % Cells& spikeph
  nunit   = max(dat.CellsInfo(:,1));
  Cells   = dat.CellsInfo(:,1);
  spikeph = zeros(nunit,3);
  spikeph(Cells,:) = dat.CellsInfo;
  % save
  dat.Cells   = Cells;
  dat.spikeph = spikeph;
  
  % SessionNP
  nlapr = numel(unique(dat.spikeRtr))-1; % first trial (0) is not a trial...
  nlapl = numel(unique(dat.spikeLtr))-1; % first trial (0) is not a trial...
  nlap = numel(unique(dat.spikeRtr)) + numel(unique(dat.spikeLtr));
  SessionNP = []; % will expand
  % for right laps
  for ilap = 1:nlapr
    % start time estimate
    currspind = dat.spikeRtr==ilap;
    currlapst = dat.spiket(currspind);
    starttime = (min(currlapst) ./ fsspike) - 0.010; % add 10ms
    % nose poke estimate, at pos = 0.1  and end time as max pos
    [dum ind] = min(abs(dat.spikeX(currspind)-0.1));
    nptime    = currlapst(ind) ./ fsspike;
    [dum ind] = max(dat.spikeX(currspind));
    endtime   = (currlapst(ind) ./ fsspike) + 0.010; % add 10ms
    SessionNP(end+1,1:4) = [starttime nptime endtime 1]; % 1 = right
  end
  % for left laps
  for ilap = 1:nlapl
    % start time estimate
    currspind = dat.spikeLtr==ilap;
    currlapst = dat.spiket(currspind);
    starttime = (min(currlapst) ./ fsspike) - 0.010; % add 10ms
    % nose poke estimate, at pos = 0.1  and end time as max pos
    [dum ind] = min(abs(dat.spikeX(currspind)-0.1));
    nptime    = currlapst(ind) ./ fsspike;
    [dum ind] = max(dat.spikeX(currspind));
    endtime   = (currlapst(ind) ./ fsspike) + 0.010; % add 10ms
    SessionNP(end+1,1:4) = [starttime nptime endtime 2]; % 2 = left
  end
  % sort by time
  [dum sortind] = sort(SessionNP(:,1));
  SessionNP = SessionNP(sortind,:);
  % save
  dat.SessionNP = SessionNP;
  
  % whlrl_ex
  whlrld = load([sessionpath 'EE.188.whl']);
  [sortedspt sortind] = sort(dat.spiket);
  sortedspx           = dat.spikeX(sortind); % get X at sorted spike times
  [whlind ind1 ind2]  = unique(round((sortedspt)./512)+1); % convert spike times to whl samples
  whlrld(whlind,7)    = sortedspx(ind1)+1; % raise position by 1 to confer to other datasets
  % set end x to 2 as this is close, but no cigar
  whlrld(round(dat.SessionNP(:,3).*fsframe),7) = 2;
  % save
  dat.whlrld = whlrld;
end
if ~isfield(dat,'whlrld')
  error(['whlrl not found for ' sessname])
end

%%%%%%%%%% Create trl
% define trials
% specify trl in time stamps at desired resolution
trl = [];
switch trialdeftype
  
  case 'mazerun'
    %%%%%%%%
    laptrl  = round(dat.SessionNP(:,[2 3]) * fsframe); % convert to frames of movie, the only time I have lindist
    nlap    = size(laptrl,1);
    lindist = dat.whlrld(:,7);
    % check for zeros at start, and fix...
    if any(lindist(laptrl(:,1))==0)
      for ilap = 1:nlap
        if lindist(laptrl(ilap,1))==0
          currbegsmp = laptrl(ilap,1);
          % find closest non-zero backwards
          pastsmp = find(lindist(1:currbegsmp),1,'last');
          futsmp  = find(lindist(currbegsmp:end),1,'first') + currbegsmp-1;
          if lindist(pastsmp)>lindist(futsmp)
            lindist(currbegsmp) = lindist(futsmp); % current position cannot be estimated, use future sample
          else
            intld   = linspace(lindist(pastsmp),lindist(futsmp),numel(pastsmp:futsmp));
            lindist(currbegsmp) = intld(currbegsmp-pastsmp+1);
          end
          smpdist = min(futsmp-currbegsmp+1,currbegsmp-pastsmp+1);
          warning(['lindist == 0 at start of lap ' num2str(ilap) ', fixing sparsely by interpolation, closest point ' num2str(smpdist) ' samples away'])
        end
      end
    end
    % correct end, might be a sample or some off of and be 0
    % if bigger, fix...
    for ilap = 1:nlap
      currendsmp = laptrl(ilap,2);
      if lindist(currendsmp) == 0
        % find closest non-zero (end of track) in lindist
        if ilap==nlap
          nonzeroind = find(lindist ~= 0);
        else
          nonzeroind = find(lindist(1:(laptrl(ilap+1,1)-1)) ~= 0); % ensure start of next lap is not used for search
        end
        [dum, ind] = min(abs(nonzeroind-currendsmp));
        newendsmp  = nonzeroind(ind);
        if abs(newendsmp-currendsmp)<5
          % just use new end, no warning
        else
          warning(['end position inaccurate in lap ' num2str(ilap) ', taking closed endpoint assumed belonging to trial, distance of ' num2str(newendsmp-currendsmp) ' samples'])
        end
        laptrl(ilap,2) = newendsmp;
      end
    end
    %     % remove laps if it doesn't contain every segment as specified below (mindist:(maxdist-mindist)/16:maxdist)
    %     remind = [];
    %     for ilap = 1:nlap
    %       beglapsmp = laptrl(ilap,1);
    %       endlapsmp = laptrl(ilap,2);
    %       lapdist   = lindist(beglapsmp:endlapsmp);
    %       mindist = 1.1;
    %       maxdist = 2;
    %       segdist   = mindist:(maxdist-mindist)/16:maxdist;
    %       % check for full lap being present
    %       if lapdist(end)<segdist(end-1) || lapdist(1)>segdist(2)
    %         % full lap not present
    %         remind = [remind ilap];
    %       end
    %     end
    %     laptrl(remind,:) = [];
    %     if ~isempty(remind)
    %       warning(['removed ' num2str(numel(remind)) '/' num2str(nlap) ' laps due to full lap not being present in lindist'])
    %     end
    nlap = size(laptrl,1);
    % fix lindist in the lap interval, which sometimes contains zeros (due to something) remove zeros by linear interpolation to get best estimate of rat position
    for ilap = 1:nlap
      beglapsmp = laptrl(ilap,1);
      endlapsmp = laptrl(ilap,2);
      currlindist = lindist(beglapsmp:endlapsmp);
      zeroind = currlindist==0;
      % find start/end of zero segments
      zsegbeg = find(diff(zeroind)==1);
      zsegend = find(diff(zeroind)==-1) + 1; % they'll always be equal size, as the end and start are always non-zero (enforced in the above laptrl fixing)
      nzseg = numel(zsegbeg);
      for izseg = 1:nzseg
        currseg = currlindist(zsegbeg(izseg):zsegend(izseg));
        currseg = linspace(currseg(1),currseg(end),numel(currseg));
        currlindist(zsegbeg(izseg):zsegend(izseg)) = currseg;
      end
      % put back fixed lindist
      lindist(beglapsmp:endlapsmp) = currlindist;
    end
    
    % loop over laps, further segment using lindist
    for ilap = 1:nlap
      % get lap start/end
      beglapsmp = laptrl(ilap,1);
      endlapsmp = laptrl(ilap,2);
      
      % create segments
      mindist = 1;
      maxdist = 2;
      segdist   = mindist:(maxdist-mindist)/(nsegments+1):maxdist;
      begsegdist = segdist(1:1:end-2)';
      endsegdist = segdist(3:1:end)';
      lapdist = lindist(beglapsmp:endlapsmp);
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
      
      % create laptrl
      begsample = begsegsmp;
      endsample = endsegsmp;
      offset    = zeros(nseg,1);
      lapnum    = ones(nseg,1) .* ilap;
      segnum    = (1:nseg).';
      trialtype = ones(nseg,1) .* dat.SessionNP(ilap,4);
      segtime   = (endsegsmp-begsegsmp) ./ fsframe;
      lapsegtrl = [begsample endsample offset lapnum segnum trialtype segtime begsegdist endsegdist];
      
      
      % add to trl
      trl(end+1:end+nseg,1:size(lapsegtrl,2)) = lapsegtrl;
    end
    
    
    
  case '5secseg'
    
    % get total length of whl file
    nsmpwhl = size(dat.whlrld,1);
    % subdivide in n sec bits
    segsmp = 1:round(5*fsframe):nsmpwhl;
    % make even
    segsmp = segsmp(1:end-mod(numel(segsmp),2));
    % create ingredients
    begsample = segsmp(1:2:end)';
    endsample = segsmp(2:2:end)';
    offset    = ones(numel(begsample),1) .* 0;
    % create trl
    trl = [begsample endsample offset];
    
    
  otherwise
    error('undefined trialdeftype')
end % switch
ntrial = size(trl,1);
%%%%%%%%%% Create trl


%%%% Debug plotting of trial definition
if plotresult
  % first plot the whole experimental time course
  figure('numbertitle','off','name','trial definition result fig 1');
  hold on
  % get all ingredients
  lindistold = dat.whlrld(:,7);
  x = [trl(:,1) trl(:,1) trl(:,2) trl(:,2)]';
  x = x ./ fsframe;
  y = trl(:,5);
  y = (y ./ max(y)) * maxdist;
  y = [zeros(ntrial,1) y y zeros(ntrial,1)]';
  % plot it all
  plot((1:numel(lindist))./fsframe,lindist);
  plot((1:numel(lindistold))./fsframe,lindistold);
  line(x,y,'color','r');
  legend('linearized distance','orig linearized distance','trials, start to end')
  xlabel('time (s)')
  ylabel('linearized distance (au)')
  
  % then plot lap specific scatters
  figure('numbertitle','off','name','trial definition result fig 2');
  nlap = numel(unique(trl(:,4)));
  for ilap = 1:nlap
    subplot(ceil(sqrt(nlap)),ceil(sqrt(nlap)),ilap)
    hold on
    % find start/end
    begsmp = round(dat.SessionNP(ilap,1) * fsframe); % start at beginning of nose poke
    endsmp = trl(find(trl(:,4)==ilap,1,'last'),2);
    % create scatter ingredients
    X = dat.whlrld(begsmp:endsmp,1);
    currlld   = lindist(begsmp:endsmp);
    Y = dat.whlrld(begsmp:endsmp,2);
    % create line ingredients
    currind = find(trl(:,4)==ilap);
    nseg    = numel(currind);
    x = [lindist(trl(currind,1)) lindist(trl(currind,1)) lindist(trl(currind,2)) lindist(trl(currind,2))]';
    y = trl(currind,5);
    y = (y ./ max(y)) * 600; % 650 is about the max of X/Y
    y = [zeros(nseg,1) y y zeros(nseg,1)]';
    % plot it all
    scatter(currlld,X,'r.')
    scatter(currlld,Y,'gr.')
    xlabel('lindist')
    ylabel('X and Y')
    line(x,y,'color','r');
    axis tight
    set(gca,'xlim',[1 2],'ylim',[0 650])
    title(['lap ' num2str(ilap)])
  end
end
%%%% Debug plotting of trial definition



%%%%%%%%%% Create spiketrains
% create trial and time cell-array containing spike trains
% select only pyramidal cells
%unitsel = sort(unique(dat.spikeind));
unitsel = sort(dat.Cells(:,1)); % use selecting by original experimenters
nunit   = numel(unitsel);
trial   = cell(1,ntrial);
time    = cell(1,ntrial);
% get timestamp samples at new sampling rate
tsssnewfs = round(dat.spiket .* spikefsfac);
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
  currtrialuid = dat.spikeind(ind);
  
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
  label{iunit} = ['unit' num2str(curruid) '_' num2str(dat.spikeph(curruid,2))];
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


% plot per trail the x-y position of the animal
figure
nlap = numel(unique(trl(:,4)));
for ilap = 1:nlap
  subplot(ceil(sqrt(nlap)),ceil(sqrt(nlap)),ilap)
  hold on
  % find start/end
  begsmp = trl(find(trl(:,4)==ilap,1),1);
  endsmp = trl(find(trl(:,4)==ilap,1,'last'),2);
  X = dat.whlrld(begsmp:endsmp,1);
  Y = dat.whlrld(begsmp:endsmp,2);
  lindist = dat.whlrld(begsmp:endsmp,7);
  scatter(lindist,X,'r.')
  scatter(lindist,Y,'gr.')
  xlabel('lindist')
  ylabel('X and Y')
  xl = dat.whlrld(round(dat.SessionNP(ilap,2) * fsframe),7);
  xl = [xl xl];
  yl = [min([X(:); Y(:)]) max([X(:); Y(:)])];
  line(xl,yl,'color','cy')
  axis tight
  set(gca,'xlim',[1 2])
end


% plot per original trail
figure
nlap = size(dat.SessionNP,1);
for ilap = 1:nlap
  subplot(ceil(sqrt(nlap)),ceil(sqrt(nlap)),ilap)
  %subplot(2,3,ilap)
  hold on
  % find start/end
  begsmp = round(dat.SessionNP(ilap,1) .* fsframe);
  endsmp = round(dat.SessionNP(ilap,3) .* fsframe);
  X = dat.whlrld(begsmp:endsmp,1);
  Y = dat.whlrld(begsmp:endsmp,2);
  lindist = dat.whlrld(begsmp:endsmp,7);
  scatter(lindist,X,'r.')
  scatter(lindist,Y,'gr.')
  xlabel('norm. position, start = SessionNP col1, end = SessionNP col3 ')
  ylabel('X and Y')
  xl = dat.whlrld(round(dat.SessionNP(ilap,2) * fsframe),7);
  xl = [xl xl];
  yl = [min([X(:); Y(:)]) max([X(:); Y(:)])];
  line(xl,yl,'color','cy')
  axis tight
  set(gca,'xlim',[1 2])
  legend('X and norm. pos.','Y and norm. pos.','end nose poke');
  title(['trial ' num2str(ilap) ' - ' sessname])
end


% plot per original trail
figure
nlap = size(dat.SessionNP,1);
for ilap = 1:nlap
  %subplot(ceil(sqrt(nlap)),ceil(sqrt(nlap)),ilap)
  subplot(ceil(sqrt(nlap)),ceil(sqrt(nlap)),ilap)
  hold on
  % find start/end
  begsmp = round(dat.SessionNP(ilap,1) .* fsframe);
  endsmp = round(dat.SessionNP(ilap,3) .* fsframe);
  X = dat.whlrld(begsmp:endsmp,1);
  Y = dat.whlrld(begsmp:endsmp,2);
  lindist = dat.whlrld(begsmp:endsmp,7);
  scatter(lindist,X,'r.')
  scatter(lindist,Y,'gr.')
  xlabel('norm. position, start = SessionNP col1, end = SessionNP col3 ')
  ylabel('X and Y')
  xl = dat.whlrld(round(dat.SessionNP(ilap,2) * fsframe),7);
  xl = [xl xl];
  yl = [min([X(:); Y(:)]) max([X(:); Y(:)])];
  line(xl,yl,'color','cy')
  axis tight
  set(gca,'xlim',[1 2])
  legend('X and norm. pos.','Y and norm. pos.','end nose poke');
  title(['trial ' num2str(ilap) ' - ' sessname])
end

% plot x and trial def of SessionNP
figure('numbertitle','off','name','trial definition result fig 1');
hold on
% get all ingredients
lindistold = dat.whlrld(:,7);
trl = round(dat.SessionNP(31,[2 3]).*fsframe);
ntrial = size(trl,1);
x = [trl(:,1) trl(:,1) trl(:,2) trl(:,2)]';
x = x ./ fsframe;
y = (1:ntrial)';
y = (y ./ max(y)) * 2;
y = [zeros(ntrial,1) y y zeros(ntrial,1)]';
% plot it all
plot((1:numel(lindistold))./fsframe,lindistold);
line(x,y,'color','r');
legend('linearized distance','trials, start to end')
xlabel('time (s)')
ylabel('linearized distance (au)')





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








% construct whlrld
fn = [sessionpath sessname '.whlrl'];
whlrld = load(fn);
figure;
hold on
for isel = 1:2
  selind = dat.whlrld(:,5) == isel;
  ind    = find(diff(selind));
  begsmp = ind(1:2:end)+1;
  endsmp = ind(2:2:end);
  for ilap = 1:numel(begsmp)
    currx = dat.whlrld(begsmp(ilap):endsmp(ilap),1);
    curry = dat.whlrld(begsmp(ilap):endsmp(ilap),2);
    currx(currx==-1) = NaN;
    curry(curry==-1) = NaN;
    currx = currx-currx(1);
    curry = curry-curry(1);
    dist = sqrt(currx.^2 + curry.^2); % distance from starting point
%     dist(isnan(dist)) = 0;
%     dist = cumsum(dist);
%     dist = dist ./ dist(end);

    dist = abs(diff(dist)); % absolute moment-by-moment change of distance from starting point
    dist(isnan(dist)) = 0;
    dist = [0; cumsum(dist)];
    dist = dist ./ dist(end);
    plot(dist)
    whlrld(begsmp(ilap):endsmp(ilap),7) = dist;
  end
end

x = whlrld(:,1);
y = whlrld(:,2);
dist = whlrld(:,7);
figure
hold on
scatter(dist,x)
scatter(dist,y,'r')

plotyy(1:numel(x),x,1:numel(x),dist)
plotyy(1:numel(x),y,1:numel(x),dist)

yfix = y;
xfix = x;

xedge = floor(min(x)):ceil(max(x));
medy1 = NaN(1,numel(xedge)-1);
medy2 = NaN(1,numel(xedge)-1);
for ix = 1:numel(xedge)-1
  %
  xind = (x > xedge(ix)) & (x < xedge(ix+1)) & y<(max(y)-min(y));
  medy1(ix) = median(y(xind));
  yfix(xind) = median(y(xind));
  %
  xind = (x > xedge(ix)) & (x < xedge(ix+1)) & y>=(max(y)-min(y));
  medy2(ix) = median(y(xind));
  yfix(xind) = median(y(xind));
end

yedge = floor(min(y)):ceil(max(y));
medx = NaN(1,numel(yedge)-1);
for iy = 1:numel(yedge)-1
  yind = (y > yedge(iy)) & (y < yedge(iy+1));
  medx(iy) = median(x(yind));
  xfix(yind) = median(x(yind));
end

  

figure
hold on
plot(dat.whlrld(:,7))
plot(whlrld(:,7)+1,'r')


figure
hold on
plot(dat.whlrld(:,7))
plot(whlrld(:,7)+1,'r')




















