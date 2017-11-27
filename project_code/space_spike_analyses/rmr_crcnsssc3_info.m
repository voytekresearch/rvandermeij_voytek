function info = rmr_crcnsssc3_info


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
info.datapath = [pathprefix 'common/data2/CRCNS_SSC3/'];
info.savepath = [pathprefix 'roevdmei/crcnsssc3/'];
if info.ontscc
  info.scratchpath = '/oasis/tscc/scratch/roevdmei/';
else
  info.scratchpath = info.savepath;
end

% sessions
info.session = {'DataSet1','DataSet2','DataSet3','DataSet4','DataSet5','DataSet6','DataSet7','DataSet8','DataSet9','DataSet10','DataSet11','DataSet12','DataSet13','DataSet14','DataSet15','DataSet16','DataSet17','DataSet18','DataSet19','DataSet20','DataSet21','DataSet22','DataSet23','DataSet24','DataSet25'};












