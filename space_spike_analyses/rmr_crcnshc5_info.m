function info = rmr_crcnshc5_info


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
info.datapath = [pathprefix 'common/data2/CRCNS_HC5/'];
info.savepath = [pathprefix 'roevdmei/crcnshc5/'];
if info.ontscc
  info.scratchpath = '/oasis/tscc/scratch/roevdmei/';
else
  info.scratchpath = info.savepath;
end

% sessions
info.session = {'i01_maze06.002','i01_maze05.005'};

