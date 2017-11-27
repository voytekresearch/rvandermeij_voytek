% 
% This script can be used to read in data from a single subject using FieldTrip.
% This script depends on the function rmr_predfaceval_definetrials.m, which should
% be on the MATLAB path. This function provides all the necessary information to
% segment the raw recordings into trials (with some control over prestim/poststim periods).
% Syncing/checking the timing provived by the EDF/BESA photo diode and PTB matfile is performed 
% in this function as well, see it for further details.
%  
% For all trials that are extracted, t=0 is the onset of the CUE, and determines the start of the trial.
% The end of each trial is determined by the onset of the FACE. The prestimulus period (prior to CUE onset) and 
% poststimulus period (after FACE onset) that are additionally extracted can be defined using the options 
% below. If both periods are 0, the trial starts at the CUE onset and ends at the FACE onset. 
%
% The output will contain a 'trialinfo' field, which contains trial-specific information, such as 
% the time from t=0 of the onset of the FACE. For more info see below.
%
% Reading in the raw recordings and segmenting it into trials is performed using FieldTrip.
% FieldTrip can be found on http://www.fieldtriptoolbox.org/
% For additional information on preprocessing options and, see the tutorials in the Preprocessing
% section of the FieldTrip tutorial list, which can be found at: http://www.fieldtriptoolbox.org/tutorial 
%
% The MATLAB-structure 'data' produced below is a FieldTrip style data-structure. 
% It will contain, among others, the following important fields:
%            hdr: [1x1 struct]   - the original header information from the BESA/EDF file
%          label: {143x1 cell}   - Nchanx1 cell-array, containing the channel names (from the BESA/EDF file)
%           time: {1x60 cell}    - 1xNtrial cell-array, containing the time axis of each trial, to be used together with data.trial
%          trial: {1x60 cell}    - 1xNtrial cell-array, containing the data of each trial, as a NchanxNsamples matrix
%        fsample: 5000           - sampling rate of the data in Hz
%     sampleinfo: [60x2 double]  - the original sample numbers of each trial on disk
%      trialinfo: [60x9 double]  - NtrialxNvariable trial-specific info variables 
%            cfg: [1x1 struct]   - the cfg used to obtain the data
%
% The trialinfo field contains the following information in each column, each row being a trial:
%  1) pred/unpred - predictive (1) or unpredictive (2) cue trial
%  2) fear/neut   - fearful (1) or neutral (2) face trial
%  3) hit/miss    - detected valence correctly (1) or incorrectly (0)
%  4) RT          - reaction time in seconds
%  5) faceonset   - onset of face stimulus in seconds (from cue onset, t=0)
%  6) respcueons  - estimated onset of response cue in seconds (from cue onset, t=0)
%  7) faceindex   - index of shown face as defined in MGL task script
%  8) cuediodur   - duration of cue on the screen as measured by the photo diode
%  9) facediodur  - duration of face+mask on the screen as measured by the photo diode
%

% set the direct path to the BESA/EDF file and PTB .mat file
datapath = '/projects/ps-voyteklab/common/data1/pred_face_val/IR35/';
datafn   = [datapath '2016021912_0003.besa'];
eventfn  = [datapath 'fp_task_19-Feb-2016_3.mat'];


% specificy options necessary obtaining segmentation details
cfg = []; % start with an empty cfg
cfg.datafile  = datafn; % for BESA: *.besa, for EDF, *.edf
cfg.eventfile = eventfn; % the EDF/BESA file and PTB mat-file should belong to the same session
cfg.diodechan = 'DC01'; % this is the label, found in the header of the EDF/BESA file, of the channel that contains 
                                % the signal from the photo diode. This is subject-specific, but is likely always one 
                                % of the first 4 channels. To obtain this name, one case use the above cfg in 
                                % ft_databrowser(cfg), selecting the first 4 channels, and observe which one contains 
                                % stepwise signals.
cfg.prestim   = 1; % the period, in seconds, before CUE ONSET that is additionally cut out (t=0 will remain CUE ONSET)
cfg.poststim  = 1; % the period, in seconds, after FACE ONSET that is additionally cut out (FACE ONSET is kept in data.trialinfo, see above)
cfg.debugflg  = true; % create a figure showing diode event detection output
trl = rmr_predfaceval_definetrials(cfg); % obtain the trl matrix, which contains the segmentation details

% read in data, using only the trl, and the BESA/EDF datafile
cfg = []; % start with an empty cfg
cfg.datafile  = datafn; % 
cfg.trl       = trl;
data = ft_preprocessing(cfg); % obtain segmented data


% After doing the above for all avialable sessions (producing a separate data-structure each),
% they can be combined by doing the following (data1, data2, etc refer to the data from each session)
datacomb = ft_appenddata([], data1, data2, data3); % in the case of 3 sessions








