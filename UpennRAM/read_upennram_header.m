function [hdr] = read_upennram_header(headerfile)

% READ_UPENNRAM_HEADER reads the header from UPenn's RAM (Restoring Active Memory) publically
% available datasets.
% 
% UPenn's RAM project is a DARPA funded project, headed by Daniel Rizzuto and Michael Kahana.
% From the RAM website: (http://memory.psych.upenn.edu/RAM)
% "The goal of RAM is to develop a fully implantable device that can electrically stimulate 
% the brain to improve memory function. The program?s immediate focus is to deliver new treatments 
% for those who have experienced a traumatic brain injury, such as veterans returning from combat. 
% In the long term, such therapies could help patients with a broad range of ailments, from Alzheimer's
% to dementia. RAM is part of a broader portfolio of programs within DARPA that support President Obama's 
% BRAIN initiative."
% Intracranial human EEG has been recorded in various tasks and is publically available at the above 
% website.
% 
% From the press release (https://news.upenn.edu/news/penns-restoring-active-memory-project-releases-extensive-human-brain-dataset):
% "The recently released dataset includes information from 700 sessions, and for every patient, 
% intracranial recording files from 100 to 200 electrode channels, neuro-anatomical information 
% indicating the location of each electrode, precise records of patient behavior and the experimental 
% design documents. To receive the raw dataset, interested researchers may request access through 
% the RAM website."
% 
% 
% Use as
% [hdr] = read_upennram_header(headerfile)
%
%      Where,
%         headerfile = string, path+name of JSON file pointing to details of a RAM dataset 
%                           IMPORTANT: this NEEDS to be the index.json file, in $PATIENTID/experiments/$TASKID/sessions/$#/behavioral/current_processed/
%                                      of the task and session of interest. This file allows to unamibigously identify the 'true' headers files:
%                                      (1) $PATIENTID/experiments/$TASKID/sessions/$#/ephys/current_processed/sources.json (binary format, nsample, etc)
%                                      (2) $PATIENTID/localizations/$#/montages/$#/neuroradiology/current_processed/contacts.json (channel names)
%                            (It is key that #2 reflects the montage used on the day of recording, incase the patients had a reimplantation, 
%                             which is ensured by using index.json.)
%      And,
%                hdr = matlab structure, containing details of the respective dataset, such as channel names, binary format, number of samples, etc
%
% 
% The read_upennram_event and read_upennram_header functions currently depends on JSONlab (by Qianqian Fang)
% from the MATLAB File Exchange, to parse the above JSON formatted files.
%
% See also READ_UPENNRAM_EVENT, READ_UPENNRAM_DATA

% Copyright (C) 2016-2017, Roemer van der Meij

% obtain sources.json info
ind       = strfind(headerfile,'behavioral');
sourcesfn = [headerfile(1:ind-1) 'ephys' headerfile(ind+10:end-10) 'sources.json'];
hdr1      = loadjson(sourcesfn);
fieldname = fieldnames(hdr1);
hdr1      = hdr1.(fieldname{:}); % failsafe in case of event that two main fields are present


% obtain patient path
ind     = strfind(headerfile,'experiments');
patpath = headerfile(1:ind-1);

% parse contents of 'hook' headerfile
hookhdr = loadjson(headerfile);

% obtain contacts.json info
contactsfn = hookhdr.info.contacts;
ind        = strfind(contactsfn,'localizations');
contactsfn = [patpath contactsfn(ind:end)];
hdr2       = loadjson(contactsfn);
fieldname  = fieldnames(hdr2);
hdr2       = hdr2.(fieldname{:}); % failsafe in case of event that two main fields are present

%%% Parse obtained information and create standardized hdr structure

% obtain list of possible channels and their data-file file extensions, combine into filenames 
label     = fieldnames(hdr2.contacts);
fileext   = cell(size(label));
for ichan = 1:numel(label)
  currind = hdr2.contacts.(label{ichan}).channel;
  fileext{ichan} = num2str(currind,'%03.0f');
end

% get number of bytes per element of dataformat (lame..., there should be matlab function for this)
tmp = cast(1,hdr1.data_format);
tmp = whos('tmp');
nbytes = tmp.bytes;

% construct header
hdr = [];
hdr.nChans       = numel(label);
hdr.nSamples     = hdr1.n_samples;
hdr.Fs           = hdr1.sample_rate;
hdr.label        = label;
hdr.labelfileext = fileext;
hdr.nSamplesPre  = 0;
hdr.nTrials      = 1;
hdr.dataformat   = hdr1.data_format;
hdr.nBytes       = nbytes;
hdr.orig1        = hdr1;
hdr.orig2        = hdr2;

    
    