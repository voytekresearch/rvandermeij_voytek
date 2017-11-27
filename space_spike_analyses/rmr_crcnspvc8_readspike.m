function data = rmr_crcnspvc8_readspike(sessionpath,fsample,trialdeftype)

% This function reads and processes the spike data from CRCNS PVC8 datasets
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
% Trialdeftype  = 'stimonly';
%
% Trialinfo fields
%   1) iimage      = image #
%   2) irep        = repetition # of image
%   3) natgrat     = natural (1) or grating (2) image
%   4) smallbig    = small (1) or big (2) natural image (NaN for grating)
%   5) nattype     = natural image category 1 to 9 (1:4 strong, 5:8 weak, 9 no dominant orientation content (3) (NaN for gratings) (within each category: 0/45/90/135 degrees)
%   6) grattype    = gratings for spat. freq. tuning (1), orientation tuning (2), size tuning (3) (NaN for natural images)
%   7) itrial      = original trial count
%
%
%
%
%

% input and error checking
if rem(fsample,1000) ~= 0
  error('fsample should be an integer multiple of 1000Hz')
end

% set paths/filenames and load data
%sessionpath = '/Users/roemer/Work/Data/tmpdata/CRCNS_PVC8/01.mat';
dat = load(sessionpath);
%%%%%%%%%% from m file accompanying dataset
%*** pre-processing: shift by 50ms to account for response latency
latency=50;
tmp = dat.resp_train;
tmp(:,:,:,1:end-latency) = dat.resp_train(:,:,:,latency+1:end);
tmp(:,:,:,end-latency+1:end) = dat.resp_train_blk(:,:,:,1:latency);
dat.resp_train = tmp;
%%%%%%%%%% from m file accompanying dataset


% parse some basics
unitind = find(dat.INDCENT);
nunit   = numel(unitind);
ntrial  = size(dat.resp_train,2) .* size(dat.resp_train,3);
ntime   = size(dat.resp_train,4);

%%%%%%%%%% Create trial matrices and trialinfo field
trial      = cell(1,ntrial);
trialinfo  = NaN(ntrial,7);
spktrain   = permute(dat.resp_train,[1 4 2 3]);
switch trialdeftype
  
  case 'stimonly'
    count = 0;
    for iimage = 1:size(spktrain,3)
      for irep = 1:size(spktrain,4)
        count = count+1;
        
        % fetch trial spiketrian
        trial{count} = spktrain(unitind,:,iimage,irep); % it will be made sparse double later
        
        % add info to trialinfo
        if iimage <= 540 % natural
          natgrat  = 1; 
          if rem(iimage,2) == 1;
            smallbig = 1;
          else
            smallbig = 2;
          end
          nattype  = ceil(rem(iimage,18)/2);
          if nattype==0
            nattype = 9;
          end
          grattype  = NaN;
        else % grating
          natgrat  = 2;
          smallbig = NaN;
          nattype  = NaN;
          if (iimage-540)<=128
            grattype  = 1;
          elseif (iimage-540)>128 &&  (iimage-540)<=(128+64)
            grattype  = 2;
          else
            grattype  = 3;
          end
        end
        trialinfo(count,:) = [iimage irep natgrat smallbig nattype grattype count]; 
      end
    end
    
  otherwise
    error('undefined trialdeftype')
end % switch
%%%%%%%%%% 


%%%%%%%%%% Create label field
label = [];
for iunit = 1:nunit
  curruid = unitind(iunit);
  label{iunit} = ['unit' num2str(curruid)];
end
%%%%%%%%%%



%%%%%%%%%% Create data structure
data            = [];
data.hdr.Fs     = fsample;
data.label      = label;
data.time       = repmat({[0 105]/1000},[1 ntrial]);
data.trial      = trial;
data.fsample    = fsample;
data.sampleinfo = [(1:ntime:(ntime*ntrial))' (ntime:ntime:(ntime*ntrial))' zeros(ntrial,1)];
data.trialinfo  = trialinfo;
data.cfg        = [];
%%%%%%%%%%



%%%%%%%%%% upsample to fsample
for itrial = 1:numel(data.trial)
  currtrial = data.trial{itrial};
  % obtain unit id's and sample id's of spikes
  [unitid, smpid] = ind2sub(size(currtrial),find(currtrial));
  smpid = smpid.*(fsample/1000);
  currtrial = sparse(unitid,smpid,ones(1,numel(smpid)),numel(data.label),20000);
  % put back
  data.trial{itrial} = currtrial;
end
data.fsample = 20000;
%%%%%%%%%%
































