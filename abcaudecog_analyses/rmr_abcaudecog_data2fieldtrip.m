function data = rmr_abcaudecog_data2fieldtrip(cfg)

%
%
% Read data from mat files into fieldtrip format (which mostly is just a parse and combining of individual channels)
% (meant to be preprocessed and segmented afterwards)
%
% Input:
%  cfg.subject    = subject name
%  cfg.session    = session name
%  PROBABLY NEED TO ADD SOME CHANNEL SELECTION, BUT WILL WAIT TILL I HAVE THE ***KEY***
%
%
%
%
%


% get details
info = rmr_abcaudecog_info;


%%%%%% Do prework
% get sampling rate as specified by ABC
switch cfg.subject
  case {'GP14','GP15','GP22','GP28','GP35','ST1','ST6','ST8'}
    fsample = 3051.76;
  case {'JH2','JH13'}
    fsample = 1000;
  otherwise
    error('woopsie')
end

% get channel file list, and don't assume each channel is always present (there could be missing channel files)
% file format JH2, JH11, JH13: gdat_XX_NAME.mat
% file format other subjkects: gdat_XX.mat
chanfiles = dir([info.datapath cfg.subject '/data/' cfg.session '/' 'gdat*']);
chanfiles = {chanfiles.name};
chanfiles = chanfiles(:);
nchan     = numel(chanfiles);
% extract channel number
channum   = NaN(size(chanfiles));
for ichan = 1:nchan
  currchan = chanfiles{ichan};
  prepos = strfind(currchan,'_');
  if ~(strcmp(cfg.subject,'JH2') || strcmp(cfg.subject,'JH11') || strcmp(cfg.subject,'JH13'))
    prepos = strfind(currchan,'_');
    postpos = strfind(currchan,'.');
  else 
    undpos = strfind(currchan,'_');
    prepos  = undpos(1);
    postpos = undpos(2);
  end
  channum(ichan) = str2num(currchan((prepos+1):(postpos-1)));
end
% sort
[dum sortind] = sort(channum);
chanfiles = chanfiles(sortind);
channum   = channum(sortind);
% report!
disp(['found ' num2str(nchan) ' channels on disk, labels ranging from ' num2str(channum(1)) ' to ' num2str(channum(end))])
%%%%%%



%%%%%% read in data and make trial
% first, read in all channels into a dat cell array
datcell = cell(nchan,1);
for ichan = 1:nchan
  fn = [info.datapath cfg.subject '/data/' cfg.session '/' chanfiles{ichan}];
  filevars = whos('-file', fn);
  if numel(filevars)>1
    error('file contains more than a single channel')
  end
  datvarname = filevars.name;
  filecontent = load(fn);
  datcell{ichan} = filecontent.(datvarname);
  clear filecontent
end
% check whether sizes are appropriate (i.e. 1 channel, same amount of channels)
datsizes = cellfun(@size,datcell,'uniformoutput',false);
datsizes = cat(1,datsizes{:});
if ~all(datsizes(:,1)==1)
  error('some channel files'' data contain more than one channel')
end
if ~all(datsizes(:,2)==datsizes(1,2))
  error('amount of samples not the same for every channel')
end
% sweet, all things are go, merge single channels and continue
dat = double(cat(1,datcell{:})); % ensure double precision
trial = [];
trial{1} = dat;
%%%%%%



%%%%%%
% correct INVERTED VOLTAGE!
% every patients ST6 and above from Stanford had an issue with sign inversion during recording, see mail from Avgusta
if strcmp(cfg.subject,'ST6') || strcmp(cfg.subject,'ST8')
    trial{1} = -1 .* trial{1};
    disp('correcting voltage inversion from Stanford recordings')
end
%%%%%%



%%%%%% make other small things
% create a 'trl'
nsamples = size(dat,2);
offset   = 0;
trl = [1 nsamples offset];

% create label 
label = [];
for ichan = 1:nchan;
  label{ichan} = ['chan' num2str(channum(ichan))];
end
label = label(:);
% create time
time = [];
time{1} = ((0:nsamples-1)-offset) ./ fsample;
%%%%%%




%%%%%% create data-matrix
data            = [];
data.hdr.Fs     = fsample;
data.hdr.nChans = numel(label);
data.hdr.label  = label;
data.label      = label;
data.time       = time;
data.trial      = trial;
data.fsample    = fsample;
data.sampleinfo = trl(:,1:2);
%data.trialinfo =  % maybe include this later
data.cfg        = [];
data.cfg.trl    = trl;
%%%%%%

















