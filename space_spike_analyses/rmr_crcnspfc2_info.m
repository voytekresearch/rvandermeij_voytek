function info = rmr_crcnspfc2_info


% get hostname for paths
[dum hostname] = system('hostname');
if ~isempty(strfind(hostname,'sdsc.edu'))
  info.ontscc = true;
else
  info.ontscc = false;
end


% paths
if info.ontscc
  pathprefix = '/projects/ps-voyteklab/';
else
  pathprefix = '/Volumes/voyteklab/';
end
info.datapath = [pathprefix 'common/data2/CRCNS_PFC2/'];
info.savepath = [pathprefix 'roevdmei/crcnspfc2/'];
if info.ontscc
  info.scratchpath = '/oasis/tscc/scratch/roevdmei/';
else
  info.scratchpath = info.savepath;
end

% sessions
info.session = {'EE.188','GG.069','EE.049','EE.152','EE.165',}; % EE.188 is used in example figs, but no behavior file...; GG.097GG./149/FF.130/FF.167/EE.088/EE.196 no spikes

