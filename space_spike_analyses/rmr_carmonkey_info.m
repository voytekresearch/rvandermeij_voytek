function info = rmr_carmonkey_info


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
info.datapath = [pathprefix 'common/data2/carmena_monkey/'];
info.savepath = [pathprefix 'roevdmei/carmonkey/'];
if info.ontscc
  info.scratchpath = '/oasis/tscc/scratch/roevdmei/';
else
  info.scratchpath = info.savepath;
end


% subj
info.subj = {'paco'};

% sessions
info.session.paco = {'paco020608c','paco020608b'};

