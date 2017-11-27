function info = rmr_upennram_info(experiments,contactsornot)


% set experiments default to FR1
if nargin == 0
  experiments   = {'FR1'}; % ,'FR2'
  contactsornot = 0;
end

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
info.datapath = [pathprefix 'common/data2/kahana_ecog_RAMphase1/'];
info.savepath = [pathprefix 'roevdmei/upennram/'];
if info.ontscc
  info.scratchpath = '/oasis/tscc/scratch/roevdmei/';
else
  info.scratchpath = info.savepath;
end

%%%
% check whether we need to use tmpdata
if ~info.ontscc && ~exist(info.datapath,'dir')
  % yep
  info.datapath = '/Users/roemervandermeij/Work/Data/tmpdata/kahana_ecog_RAMphase1/';
  warning('using temporary data')
end
%%%

%%%
% load r1.json
fnsuff = 'session_data/experiment_data/protocols/r1.json';
if exist([info.datapath fnsuff],'file')
  jsonpath = [info.datapath fnsuff];
else
  try
    jsonpath = ['/Users/roemervandermeij/Work/Data/tmpdata/kahana_ecog_RAMphase1/' 'session_data/experiment_data/protocols/r1.json']; % temp path
  catch
    error('r1.json not found')
  end
end
jsoninfo = loadjson(jsonpath);
jsoninfo = jsoninfo.protocols.r1.subjects;
% remove subjects without (cat)FRs
allsubj  = fieldnames(jsoninfo);
fieldrem = [];
for isubj = 1:numel(allsubj)
  currsubj = allsubj{isubj};
  exp      = fieldnames(jsoninfo.(currsubj).experiments);
  if ~any(ismember(exp,experiments))
    fieldrem{end+1} = currsubj;
  end
end
jsoninfo = rmfield(jsoninfo,fieldrem);
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subject names
info.subj         = fieldnames(jsoninfo);
info.subjsel1     = {'R1032D','R1128E','R1156D','R1149N','R1034D','R1162N','R1033D','R1167M','R1175N','R1154D','R1159P','R1080E','R1142N','R1059J','R1135E','R1147P','R1020J','R1045E','R1075J','R1084T','R1120E','R1129D','R1151E','R1155D','R1166D'}; % only 1000Hz R1068J (39) no frontal)
info.agesel1      = [ 19,      26,      27,      28,      29,      30,      31,      33,      34,      36,      42,      43,      43,      44,      47,      47,      48,      51,      50,      25,      33,      34,      36,      37,      38     ]; % only 1000Hz
info.subjsel2     = {'R1032D','R1006P','R1086M','R1177M','R1128E','R1156D','R1039M','R1149N','R1034D','R1112M','R1162N','R1033D','R1167M','R1102P','R1121M','R1175N','R1060M','R1089P','R1154D','R1003P','R1053M','R1066P','R1127P','R1159P','R1080E','R1142N','R1059J','R1067P','R1018P','R1135E','R1147P','R1001P','R1020J','R1002P','R1036M','R1045E','R1075J','R1084T','R1120E','R1129D','R1151E','R1155D','R1166D'}; % inc 500Hz
info.agesel2      = [ 19,      20,      20,      23,      26,      27,      28,      28,      29,      29,      30,      31,      33,      34,      34,      34,      36,      36,      36,      39,      39,      39,      40,      42,      43,      43,      44,      45,      47,      47,      47,      48,      48,      49,      49,      51,      50,      25,      33,      34,      36,      37,      38     ]; % inc 500Hz
info.subjsel3     = {'R1032D','R1128E','R1034D','R1167M','R1142N','R1059J','R1020J','R1045E'}; % R1059J good as example for ECoG cleaning
info.subjselmain  = {'R1020J','R1032D','R1045E','R1075J','R1080E','R1120E','R1128E','R1147P','R1151E','R1154D','R1166D'}; 
info.subjselmains = {'R1075J','R1020J','R1045E','R1080E','R1154D','R1032D','R1120E','R1128E','R1151E','R1147P','R1166D'}; % R1154D/R1166D extra buzz!
info.subjselexpa  = {'R1034D','R1059J','R1135E','R1142N','R1149N','R1162N'};

% extract info from jsoninfo
% format: info.$subject.(cell per locmont).(struct per session).(session paths/etc)
for isubj = 1:numel(info.subj)
  currsubj = info.subj{isubj};
  % fetch session info
  exp       = fieldnames(jsoninfo.(currsubj).experiments);
  expcount  = 0;
  session = [];
  for iexp = 1:numel(exp)
    switch exp{iexp}
      case experiments % {'catFR1','catFR2','FR1','FR2','PAL1','PAL2','YC1','YC2'} 
        jsonsessions = jsoninfo.(currsubj).experiments.(exp{iexp}).sessions;
        jsonsessname = fieldnames(jsonsessions);
        for isess = 1:numel(jsonsessname)
          expcount = expcount+1;
          session(expcount).headerfile   = ['session_data/experiment_data/' jsonsessions.(jsonsessname{isess}).task_events(1:end-16) 'index.json'];
          session(expcount).datadir      = ['session_data/experiment_data/' jsonsessions.(jsonsessname{isess}).task_events(1:end-45) 'ephys/current_processed/noreref/'];
          session(expcount).eventfile    = ['session_data/experiment_data/' jsonsessions.(jsonsessname{isess}).task_events];
          session(expcount).locmont      = [jsonsessions.(jsonsessname{isess}).localization '-' jsonsessions.(jsonsessname{isess}).montage];
          session(expcount).experiment   = exp{iexp};
        end
    end
  end
  % parse locmont
  locmont = {session.locmont};
  [unilocmont dum uniind] = unique(locmont);
  for iuniloc = 1:numel(unilocmont)
    info.sessions.(currsubj){iuniloc} = session(uniind==iuniloc);
  end
  % fetch locmont
  if contactsornot
    for iuniloc = 1:numel(unilocmont)
      contactsjson = [info.datapath 'session_data/experiment_data/protocols/r1/subjects/' currsubj '/localizations/' unilocmont{iuniloc}(1) '/montages/' unilocmont{iuniloc}(3) '/neuroradiology/current_processed/contacts.json'];
      if exist(contactsjson,'file')
        contacts = loadjson(contactsjson);
        contacts = contacts.(currsubj).contacts;
        info.montages.(currsubj){iuniloc} = contacts;
      else
        info.montages.(currsubj){iuniloc} = NaN;
        warning(['contacts.json not found for ' currsubj ': ' contactsjson])
      end
    end
  end
end
% remove specific ones with problems, lame...
%info.R1156D.session(4) = [];
%info.R1034D.session(1) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandstop frequencies found using rmr_findlinespectra (to be removed with a but, 2 order, twopass) (assuming a hardcut at 250) ALL DEPTHS ARE REMOVED! CLEANING ONLY PERFORMED ON SURFACE ELECTRODES
info.R1020J.bsfilt.peak      = [60  120 180 219.9 240 300];
info.R1020J.bsfilt.halfbandw = [0.5 0.5 0.6 0.5   0.5 0.6];
info.R1032D.bsfilt.peak      = [60  120 180 240 300.1];
info.R1032D.bsfilt.halfbandw = [0.5 0.5 0.5 0.5 0.5];
info.R1033D.bsfilt.peak      = [];
info.R1033D.bsfilt.halfbandw = [];
info.R1034D.bsfilt.peak      = [60  61.1 120 172.3 180 183.5 200 240 281.1 296.3 300 305.7];
info.R1034D.bsfilt.halfbandw = [0.5 0.5  0.5 0.5   0.5 0.5   0.5 0.6 0.5   0.8   0.6 0.5];
info.R1045E.bsfilt.peak      = []; % nothing needed
info.R1045E.bsfilt.halfbandw = []; % nothing needed
info.R1059J.bsfilt.peak      = [60  180 220 240 300];
info.R1059J.bsfilt.halfbandw = [0.5 0.5 0.5 0.5 0.5];
info.R1068J.bsfilt.peak      = [];
info.R1068J.bsfilt.halfbandw = [];
info.R1080E.bsfilt.peak      = [59.9 119.9 179.8 239.7 299.7];
info.R1080E.bsfilt.halfbandw = [0.5  0.5   0.5   0.5   0.5];
info.R1128E.bsfilt.peak      = [59.9 179.8 299.7];
info.R1128E.bsfilt.halfbandw = [0.5  0.5   0.5  ];
info.R1135E.bsfilt.peak      = [59.9 119.9 179.8 239.7 299.7];
info.R1135E.bsfilt.halfbandw = [0.5  0.5   0.5   0.5   0.5];
info.R1142N.bsfilt.peak      = [60  120 180 240 300];
info.R1142N.bsfilt.halfbandw = [0.5 0.5 0.5 0.5 0.5];
info.R1147P.bsfilt.peak      = [60 83.2 100 120 140 160 166.4 180 200 221.4 240 260 280 300];
info.R1147P.bsfilt.halfbandw = [0.5 0.5 0.5 0.5 0.5 0.5 0.5   0.5 0.5 3.6   0.5 0.7 0.5 0.5];
info.R1149N.bsfilt.peak      = [60  120 136 180 196.5 211.7 220 226.8 240 241.9 257.1 272.1 280 287.3 300];
info.R1149N.bsfilt.halfbandw = [0.6 0.5 0.5 0.8 0.5   0.5   0.5 0.5   1   0.5   0.5   0.5   0.5 0.5   0.9];
info.R1154D.bsfilt.peak      = [60  100 120 138.6 160 172.3 180 200 205.9 218.5 220 222.9 225.1 240 260 280 300];
info.R1154D.bsfilt.halfbandw = [0.5 0.5 0.5 1     0.5 0.5   0.5 0.5 0.5   0.5   0.7 2.4   0.5   0.5 0.5 0.6 0.5];
info.R1156D.bsfilt.peak      = [];
info.R1156D.bsfilt.halfbandw = [];
info.R1159P.bsfilt.peak      = [60  100.1 113.9 120 124.4 139.8 160.1 180 199.9 220 222.7 224.8 240 259.9 280.1 299.9];
info.R1159P.bsfilt.halfbandw = [0.5 0.5   0.5   0.5 0.5   0.8   0.5   0.5 0.5   1.1 0.7   0.7   0.5 0.6   0.5   0.6];
info.R1162N.bsfilt.peak      = [60  120 180 238.6 239.9 300];
info.R1162N.bsfilt.halfbandw = [0.5 0.5 0.5 0.5   0.8   0.6];
info.R1167M.bsfilt.peak      = [60  96.5 100.2 120 140.1 160 180 184.1 200 220.5 240 259.8 280.1 300];  
info.R1167M.bsfilt.halfbandw = [0.5 3    1   0.5 1.5   1   0.5 1.5   1   3     0.5 1.6   0.5   0.5];
info.R1175N.bsfilt.peak      = [60  120 180 240 280 300.1]; % [60.1 120 122.7 159.9 180.4 220  240.2 244.7 245.7 259.9 260.2 280  294.4 296 301.3];
info.R1175N.bsfilt.halfbandw = [0.5 0.5 0.8 0.8 0.5 1.5];   % [1.7  5.1 0.5   0.75  5.95  0.75 7.85  0.5   1.2   0.75  0.75  0.95 0.75  0.5 10.45];
%
info.R1075J.bsfilt.peak      = [60  120 180 220 240 300]; % [60 120 179.9 220.1 240 300 132 140 156 160 204 228 252 260 280]; (when including all L channels)
info.R1075J.bsfilt.halfbandw = [0.5 1.5 2   0.5 2   0.5]; % [1  6   8     0.5   7   1   0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]; (when including all L channels)info.R1084T.bsfilt.peak      = [];
info.R1084T.bsfilt.halfbandw = [];
info.R1120E.bsfilt.peak      = [60  179.8 299.7];
info.R1120E.bsfilt.halfbandw = [0.5 0.5   0.5];
info.R1129D.bsfilt.peak      = [];
info.R1129D.bsfilt.halfbandw = [];
info.R1151E.bsfilt.peak      = [60  100 120 123.7 140 160 180 200 210.2 215 220.1 240 247.3 259.9 300];
info.R1151E.bsfilt.halfbandw = [0.5 0.5 0.5 0.5   0.5 0.5 0.5 0.5 0.5   0.5 0.5   0.5 0.5   0.5   0.5];
info.R1155D.bsfilt.peak      = [];
info.R1155D.bsfilt.halfbandw = [];
info.R1166D.bsfilt.peak      = [60  100.1 120 140 160 180 200 218.3 220.1 223.8 240 260 280 300];
info.R1166D.bsfilt.halfbandw = [0.5 0.5   0.5 0.5 0.5 0.5 0.5 1.2   0.7   1.7   0.5 0.5 0.5 0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   


%[60  120 180 240 280 300.1]
%[0.5 0.5 0.8 0.8 0.5 1.5]

    %edgeartlen: 3.1840

     

%%%%%%%%%%%%%%%%%%%%%
% Length of edge artifact longest filter above ALL DEPTHS ARE REMOVED! CLEANING ONLY PERFORMED ON SURFACE ELECTRODES
info.R1020J.bsfilt.edgeartlen = 3.184;
info.R1032D.bsfilt.edgeartlen = 3.1050;
info.R1033D.bsfilt.edgeartlen = [];
info.R1034D.bsfilt.edgeartlen = 3.2237;
info.R1045E.bsfilt.edgeartlen = 0.2; % low pass only
info.R1059J.bsfilt.edgeartlen = 3.1840;
info.R1068J.bsfilt.edgeartlen = [];
info.R1080E.bsfilt.edgeartlen = 3.1852;
info.R1128E.bsfilt.edgeartlen = 3.1852;
info.R1135E.bsfilt.edgeartlen = 3.1852;
info.R1142N.bsfilt.edgeartlen = 3.1840;
info.R1147P.bsfilt.edgeartlen = 3.1840;
info.R1149N.bsfilt.edgeartlen = 3.0980;
info.R1154D.bsfilt.edgeartlen = 3.1840;
info.R1156D.bsfilt.edgeartlen = [];
info.R1159P.bsfilt.edgeartlen = 3.1840;
info.R1162N.bsfilt.edgeartlen = 3.1840;
info.R1167M.bsfilt.edgeartlen = 3.1840; % probably safe, had to guess certain broad filters
info.R1175N.bsfilt.edgeartlen = 3.0940;
%
info.R1075J.bsfilt.edgeartlen = 3.1840; % 3.0820;  (when including all L channels)
info.R1084T.bsfilt.edgeartlen = [];
info.R1120E.bsfilt.edgeartlen = 3.1852;
info.R1129D.bsfilt.edgeartlen = [];
info.R1151E.bsfilt.edgeartlen = 3.1840;
info.R1155D.bsfilt.edgeartlen = [];
info.R1166D.bsfilt.edgeartlen = 3.1840;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bad channels - ALL DEPTHS ARE REMOVED! CLEANING ONLY PERFORMED ON SURFACE ELECTRODES
info.R1020J.badchan     = {'RAH7','RSTB5','RSTB8','RPH7','RFB1','RFB8','RFB4','RAM8',... % broken (last few rather weakly)
                           'RAT1','RAT2','RAT3','RAT4','RAT5','RAT6','RAT7','RAT8',... % RAT1-8 spiklety, spikey
                           'RPTB1','RPTB2','RPTB3','RPTB4','RPTB5','RPTB6','RPTB7','RPTB8','RPTA1','RPTA2','RPTA3','RPTA4','RPTA5','RPTA6','RPTA7','RPTA8','RFE1','RFE2','RFE3','RFE4','RFE5','RFE6','RFE7','RFE8',... % spikey and wonky
                           'RFB2','RSTB6','RSTB7','RSTC1','RSTC2','RSTC3','RSTC4','RSTC5','RSTC6','RSTA4','RBST4','RFA8',... % spikey, buzzy, etc
                           'RAH*','RPH*'}; % depths
info.R1032D.badchan     = {'LFS8','LTS8','RTS8','LID12','LOFD12','LOTD12','RID12','ROFD12','ROTD12',... % broken
                           'LID*','LOFD*','LOTD*','RID*','ROFD*','ROTD*'}; % depths
info.R1033D.badchan     = {};
info.R1034D.badchan     = {'LFG1','LFG16','LFG24','LFG32','LFG8','LIHG16','LIHG24','LIHG8','LOFG12','LOFG6','LTS8','RIHG16','RIHG8','LOTD12','LOTD7',... % broken
                           'LFG11','LFG10','LFG12','LFG14','LIHG17','LIHG18','LOFG10','LTS4','LTS5','LTS6','LOFG8',... % some of these show huge depolarizations, very odd
                           'LOTD*'}; % depths
info.R1045E.badchan     = {'LIFS10','RPTS7','LPHD9','RPHD1','RPHD7'... % broken
                           'RAFS7',... % spikey
                           'LAHD*','LAMYD*','LMHD*','LPHD*','LPHGD*','LSTGID*','RAHD*','RAMYD*','RPHD*','RPHGD*','RSTGID*'}; % depths
info.R1059J.badchan     = {'LDC1','LDC2','LDC3','LDC4','LDC5','LDC6','LDC7','LDC8','RIHA1','LDA2','LFC1','BDC7','LDC1','LIHA1','RIHB1','LFB3','RAT1','RAT8','RDC7',... % broken
                           'LDA1','LDA2','LDA3','LAT1','LAT2','LAT3','LAT4','LAT5','LAT6','LAT7','LAT8','RAT2','RAT3','RAT6','RAT7','RPT1','LSTA1','LSTA2','LSTA3','LSTA4','LSTA5','LSTA6','LSTA7','LSTA8','LSTB1','LSTB2','LSTB3','LSTB4','LSTB5','LSTB6','LSTB7','LSTB8','LFB7','LFD1',... % spikey, buzzy (LFA4/5/6/LFB6/7 might be salvagable, borderline)
                           'LDA*','LDB*','LDC*','RDA*','RDB*','RDC*'}; % depths
info.R1068J.badchan     = {};
info.R1080E.badchan     = {'RLFS7','L9D7','R10D1','R12D7','R2D2','L11D8','R8D1','RSFS4','L5D10','R10D7',... % broken
                           'RLFS4','RPTS7',... % buzzy by itself/spikey
                          'L11D*','L13D*','L1D*','L3D*','L5D*','L7D*','L9D*','R10D*','R12D*','R14D*','R2D*','R4D*','R6D*','R8D*'};
info.R1128E.badchan     = {'RTRIGD10','RPHCD9',... % broken
                           'RANTTS1','RANTTS2','RANTTS3','RINFFS1',...  % spikey  (initially also removed RANTTS4, but not really necessary
                           'LAHCD*','LAMYD*','RAHCD*','RAMYD*','RENTD*','RFUSFD*','RINSD*','RITGD*','RLATTD*','RMTGD*','RPHCD*','RSTGD*','RTEMPD*','RTRIGD*'}; % depths
info.R1135E.badchan     = {'LHCD9','RPHCD1','RPHCD9','RSUPPS1','RSUPPS3','RSUPPS5','RSUPPS7','RSUPPS9',...  % broken  RSUPPS1/3/5/7/9 broken in sess3
                           'RANTTS3','RANTTS5','RPOSTS3','RSUPFS6','RLATPS1',...
                           'LAMYD*','LHCD*','LROI1D*','LROI2D*','LROI3D*','RAHCD*','RAMYD*','RFINSD*','RFROID*','RINS1D*','RINS2D*','RLATFD*','RPHCD*','RROI1D*','RROI2D*','RROI3D*','RROI4D*','RROI5D*'}; % depths
info.R1142N.badchan     = {'ALT6','RMIH1','RMIH2',... % broken
                           'PST1','PLT1','PLT2','PLT3','PLT4','PLT6'... % remove PST1 or not... spikelets  ,'ALT1','ALT2','ALT3','ALT4','AST1','AST2','AST3','AST4','LPIH1','LPIH2','LPIH3','LPIH4','LMIH1','LMIH2','LMIH3','LMIH4','LMIH5','LMIH6','RPPI1','LAIH1','LAIH2','LAIH3','LAIH4'
                           'AD*','PD*'}; % depths
info.R1147P.badchan     = {'LGR1','LGR29','LGR64','LDH1','LPST6','LGR*','LPT*','LSP*',... % removing LGR/LPT/LSP allows extreme noise to be referenced out
                           'LMST3','LMST4','LAST2','LAST3','LPST1',... % spikey
                           'LDA*','LDH*'}; % depth's
info.R1149N.badchan     = {'ALEX1','ALEX8','AST2',... % broken
                           'ALEX*','TT4','TT5','TT6','G16','G10','G8',...  % there's something wrong with the ALEX's
                           ''}; % no depths recorded
info.R1154D.badchan     = {'LTCG23','LTCG*',... % broken, LTCG needs to be removed in order to make line noise filterable (non-equal reference)
                           'LTCG25','LSTG8','LSTG2',... % buzzy (LSTG2/8 were to wonky)
                           'LOTD*'}; % depths
info.R1156D.badchan     = {};
info.R1159P.badchan     = {'LG1','LG16','LG24','LG25','LG30','LG33','LG34','LG35','LG36','LG38','LG48','LG56','LG62','LO5','LG49','LG64',... % broken
                           'LG3','LG4','LG5','LG6','LG7','LG8','LG9','LG1*','LG2*','LG30','LG31','LG32',...  % LG1-32 have different noise  ( 'LG*','LG3','LG4','LG5','LG6','LG7','LG8','LG9','LG1*','LG2*','LG30','LG31','LG32')
                           'LAST*','LFP*','LMST*','LOF*','LO*','LPST*','LP*','RAT*','RFP*','RIH*','RMT*','ROF*','RPT*','RP*',... % surface sets that have different noise than LG33-64, and need to be removed to be able to use LG33-64
                           'LDA*','LDH*','RDA*','RDH*'}; % depths
info.R1162N.badchan     = {'AST2',... % broken
                           'PST1','PST2','PST3','ATT4','ATT5','ILT1',... % buzz/wonkey
                           'PD*','G*','AD*'}; %depths
info.R1167M.badchan     = {'LP7','LP8','LPT19','LPT12','LP5','LAI1','LAT12','LP3','LPT6','LPT20',... % broken
                           'LAT1','LAT2','LAT3','LAT4','LAT9','LAT10','LAT11','LPT11','LPT21','LPT22','LPT23','LPT15','LPT16','LPT10','LPT9','LPT17','LPT14','LPT13','LPT7','LOF6','LOF7','LOF8','LPT8',...  % spikey
                           'LAD*','LPD*'}; % depths
info.R1175N.badchan     = {'RPST2','RPST3','RPST4','RPT6','RAT8', 'RSM6','RAF4','RAF1','RAF2',... % broken 
                           'LAT3','RAT1','RAT2','RAT3','RAT4','RAT5','RAT6','LPST1','RMF3','RAST1', 'RAST2', 'RAST3', 'RAST4','RMP2','RMP3','RMP4',... % buzzy/spikey
                           'RMF1','RMF2','RMF3','RMF4','RMF5','RPP*'}; % didnt share the main noise source ('RMP2','RMP3','RMP4','RMP5','RPT3','RPT4' ) 
info.R1075J.badchan     = {'RFD1','RFD8','LFB1','LFC1','LFD1','RFD3','RFD4','L*',... % broken. removed entire Left section, very broad noise on all electrodes, not able to reference out. Keeping them would mean removing at least ~30Hz from the 100 to 200 range, and yet still not lead to double channels, as more L channels seem dead in addition
                           'RFA8',... % buzzy/spikelets
                           ''}; % no depths
info.R1084T.badchan     = {};
info.R1120E.badchan     = {'RAMYD6','LANTS10','LPHD9','LANTS5',... % broken LANTS5 only broken in large chunk of sess2 (if it needs to be included, than the last third of sess2 needs to be removed, and the first 2/3 should be rechecked, there were about 10s that needs to be removed there (it buzzed, sortof)
                           'LANTS2','LANTS3','LANTS8','LPOSTS1',...  % spikey/wonky/spikelets I want to remove LANTS2/3/4/5, they act weird (small spikes consistent on slower wave), but, sparing them as their temporal
                           'LAHD*','LAMYD*','LFARPD*','LMHD*','LPHD*','RAHD*','RAMYD*','RFARPD*','RMHD*','RPHD*'}; % depths
info.R1129D.badchan     = {};
info.R1151E.badchan     = {'RPHD8','LOFMID1','LOFLD9','LMHD8',... % broken
                           'LMTS1','LATS10',... % spikey (maybe LMSTS1/LATS10 can be rescued if removing many segments)
                           'LAHD*','LAMYD*','LMHD*','LOFLD*','LOFMED*','LOFMID*','LPHD*','RAHD*','RAMYD*','RMHD*','ROFLD*','ROFMED*','ROFMID*','RPHD*'}; % depths
info.R1155D.badchan     = {};
info.R1166D.badchan     = {'LFPG14','LFPG15','LFPG16',...  % broken (Tammy also had LFPG10)
                           'LFPG6','LFPG7','LSFPG*',... % out of wack; removing LSFPG because of grid specific reference
                           'LID*','LPD*','LSD*'}; % depths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% info.R1175N.badchan.broken    = {'RAT8', 'RPST2', 'RPST3', 'RPST4', 'RPT6', 'RSM6', 'RAF4'};
% info.R1175N.badchan.epileptic = {'RAT2', 'RAT3', 'RAT4', 'RAT5', 'RAT6' ... % synchronous spike on bump
%     'RAT1', 'RMF3', ... % IEDs isolated to single channels
%     'LPST1', 'RAST1', 'RAST2', 'RAST3', 'RAST4' ... % more IEDs
%     };


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Presence of sample-wise artifact specification - ALL DEPTHS ARE REMOVED! CLEANING ONLY PERFORMED ON SURFACE ELECTRODES
info.R1020J.trlartfctflg = 1; % (main) later remembered vs forgotten: 107:180 trials (age 48, 37surf)
info.R1032D.trlartfctflg = 1; % (main) later remembered vs forgotten:  92:196 trials (age 19, 29surf)
info.R1033D.trlartfctflg = 0; 
info.R1034D.trlartfctflg = 1; % (expa) later remembered vs forgotten:  47:471 trials (age 29, 77surf) (very slow ECoG) 
info.R1045E.trlartfctflg = 1; % (main) later remembered vs forgotten:  81:168 trials (age 51, 45surf)  
info.R1059J.trlartfctflg = 1; % (expa) later remembered vs forgotten:  35:387 trials (age 44, 98surf)  
info.R1068J.trlartfctflg = 0; 
info.R1080E.trlartfctflg = 1; % (main) later remembered vs forgotten: 106:272 trials (age 43, 17surf)
info.R1128E.trlartfctflg = 1; % (main) later remembered vs forgotten: 140:154 trials (age 19, 19surf) 
info.R1135E.trlartfctflg = 1; % (expa) later remembered vs forgotten:  79:775 trials (age 47, 25surf)
info.R1142N.trlartfctflg = 1; % (expa) later remembered vs forgotten:  46:240 trials (age 26, 78surf)  
info.R1147P.trlartfctflg = 1; % (main) later remembered vs forgotten:  96:424 trials (age 47, 26surf)
info.R1149N.trlartfctflg = 1; % (expa) later remembered vs forgotten:  51:179 trials (age 28, 57surf)
info.R1154D.trlartfctflg = 1; % (main) later remembered vs forgotten: 245:592 trials (age 36, 32surf) 
info.R1156D.trlartfctflg = 0; 
info.R1159P.trlartfctflg = 1; % (removed: too few frontal) later remembered vs forgotten:  40:127 trials (age 42, 22surf)
info.R1162N.trlartfctflg = 1; % (expa) later remembered vs forgotten:  76:208 trials (age 30, 39surf)
info.R1167M.trlartfctflg = 1; % (removed: flat slope) later remembered vs forgotten: 162:200 trials (age 33, 35surf) (flat slope)
info.R1175N.trlartfctflg = 0; % (removed: broad ref noise after rereferencing, ton of electric blips, probably causing the former) later remembered vs forgotten: 58:220 trials (age 34, 65surf)
%
info.R1075J.trlartfctflg = 1; % (main) later remembered vs forgotten: 141:433 trials (age 50, 57surf)
info.R1084T.trlartfctflg = 0;
info.R1120E.trlartfctflg = 1; % (main) later remembered vs forgotten: 206:387 trials (age 33, 32surf)
info.R1129D.trlartfctflg = 0;
info.R1151E.trlartfctflg = 1; % (main) later remembered vs forgotten: 208:548 trials (age 36, 18surf)
info.R1155D.trlartfctflg = 0;
info.R1166D.trlartfctflg = 1; % (main) later remembered vs forgotten: 123:723 trials (age 36, 27surf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes on artifacts  - ALL DEPTHS ARE REMOVED! CLEANING ONLY PERFORMED ON SURFACE ELECTRODES
%%% Main list
% R1020J 4- Mostly buzz, the occasional big event, quite some channels removed that were spikey/wonky (quite liberal). The occasional hint of mechanical artifacts
% R1032D 3- LFS1/2 often showed large slow depolarizations, otherwise decent, not much removed
% R1045E 4- some broadband buzzing quite often, but that's about it
% R1075J 5- sess1/2: removed some minor events/buzz (sess2 had much more buzz on 5ish channels). Buzz with oscillating power spectra was hard to remove channel wise, so had use broad filters.
% R1080E 4- sess1: no major event or anything, just occasional strong buzz, some spikelets. Sess2 similar, but buzz was often focused on a single channel, occasionally big
% R1120E 3- sess1: removed some bigger events, not many, sess2: same but more, a little sortof buzz too (for line noise, last bit of each recording had to be remove due to huge step function)
% R1128E 3- RSUPPS1 often wildy oscillating, some extreme depolarizations, RINFPS3 sometimes joins in. Threw out some of the more extreme cases
% R1147P 2- sess1/2/3: removed lots of small spike-waves and buzz (sess3 had less spike waves and smaller buzz)
% R1151E 3- sess1: not really much, pretty clean, some spikey some wonky, but kept most in. Sess2: similar, removed couple of bigger spikes         (sess3 is broken after ~1490s, all trials occur before that)
% R1154D 4**- sess1: removed mostly buzz, an occasional biggish event. Sess2/3 lots more buzz.         (LTCG needed to be removed in order to make line noise filterable (non-equal reference)) (last ~100s of sess2 and first ~260s of sess3 were 'unrecorded')
% R1166D 2**- sess1: buzz buzz buzz buzz. Sess2 much less buzz. Sess3 more buzz than sess2 again. All of them also weak buzz.      (LSFPG need to be removed due to specific reference (or other) noise
% R1167M (REMOVED: flat slope not cleanable, on all channels) tricky cleaning wrt to line/broadband noise and channels. Sess1 some bigger events, some buzz removed. Small high freq noise often present, not removable... expect flatter slope. Sess2 same as 1, except for the small high freq noise.
%%% 'Expansion' list
% R1034D - (<10%) sess1: very few events removed. sess2: same, few more bigger ones removed, but session had many more synced low amplitude events. Sess3 same as sess2, less bigger ones removed. This subjects has something odd about him/her, like everything is stretched out, much much more <delta activity than normal, medication maybe? Massive depolarizations on some channels
% R1059J - (<10%) sess1: removed quite some spikes (corr: removed extra channels to save trials), the occasional buzz. Sess2 same as sess1, except larger buzz segments
% R1135E - (<10%) sess1/2/3/4: horrible cleaning, everything looked suspicious (and therefore not suspicious), many big waves/chaos moments, often together with spikey activity. Sess3, different, much less slow activity, different meds? Striking difference
% R1142N - (I thought many spikelets, but seems more mechanical and only in the beginning). Very few events thrown out during experiment, weird spikelets aggresivley thrown out
% R1149N - terribly awful cleaning, probably very inconsistent when going through again. many large synced events. Threw out many due to uncertainty, kept a lot in due to uncertainty.
% R1159P (REMOVED: too few frontal) - continuously restless, very little quietude. Removed only a handful of things.
% R1162N - removed a few bigger synced events
% R1175N - removed a few bigger events that were synced up, a f***ing ton of electric noise (jumps/buzz/etc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% OLD NOTES %%%%%%%%%%%%%%%%
% R1020J: ONLY FOR FR1 sessions (first 1 of 5). Mostly buzz, the occasional big event, quite some channels removed that were spikey/wonky (quite liberal). The occasional hint of mechanical artifacts
% R1032D: LFS1/2 often showed large slow depolarizations, otherwise decent, not much removed
% R1033D: 
% R1034D: sess1: very few events removed. sess2: same, few more bigger ones removed, but session had many more synced low amplitude events. Sess3 same as sess2, less bigger ones removed. This subjects has something odd about him/her, like everything is stretched out, much much more <delta activity than normal, medication maybe? Massive depolarizations on some channels
% R1045E: some broadband buzzing quite often, but that's about it
% R1059J: sess1: removed quite some spikes (corr: removed extra channels to save trials), the occasional buzz. Sess2 same as sess1, except larger buzz segments
% R1068J: 
% R1080E: 
% R1128E: RSUPPS1 often wildy oscillating, some extreme depolarizations, RINFPS3 sometimes joins in. Threw out some of the more extreme cases
% R1135E: 
% R1142N: (I thought many spikelets, but seems more mechanical and only in the beginning). Very few events thrown out during experiment, weird spikelets aggresivley thrown out
% R1147P: 
% R1149N: 
% R1154D: 
% R1156D: 
% R1159P: 
% R1162N: 
% R1167M: tricky cleaning wrt to line/broadband noise and channels. Sess1 some bigger events, some buzz removed. Small high freq noise often present, not removable... expect flatter slope. Sess2 same as 1, except for the small high freq noise.
% R1175N: 
%
% R1075J: Sess1/2: removed some minor events/buzz (sess2 had much more buzz on 5ish channels). Buzz with oscillating power spectra was hard to remove channel wise, so had use broad filters.
% R1084T:
% R1120E:
% R1129D:
% R1151E:
% R1155D:
% R1166D:
%%% Pre-selection notes
% **R1032D (19) - 1sess - good, cleaned
% **R1128E (26) - 1sess - medium-good, cleaned
% R1156D (27) - 3sess (1sess error) - bad - all sessions: either utterly broken or grid/strip specific reference (seems likely), hard to see without going deeper 
% R1149N (28) - 1sess - bad - all channels often sync up in spikes and wonky oscillatory patterns (feels like the normal is crazy, also the occasional hint of ref noise)
% -**R1034D (29) - 3sess - medium - sess1/2/3: decent, some spiking but not much. sess1: unfilterable continuous spike buzz (CORRECTION, incorrectly observed, removable by low-pass). Something weird about this patient, many channels too flat
% R1162N (30) - 1sess - bad-medium - lot of channels showing simultaneous wonky oscillations, quite a few big IEDs
% R1033D (31) - 4sess (FR1/2) - bad - sess1/2/3/4: many big events of all sorts, and strongly sharp wonky oscillations, 3 becomes waterslide batshit crazy, 4 continues on that path, serious electrical problems
% *R1167M (33) - 2sess - medium - sess1/2: frequent major events on some channels but somewhat isolated, some wonky oscillations but not involving many channels, some spiklets in sess2
% R1175N (34) - 1sess - medium-bad - somewhat okay, but many wonky oscillations all over. Not very strong, but not very mild either... 
% R1154D (36) - 3sess - bad - sess1/2: only one temporal grid remaining after discarding due to broadband noise, then many global wonky oscillation events. Sess3 WTF? Wasn't even recorded?
% R1068J (39) - 3sess - medium-bad - sess1/2/3: little odd, but probably salvagable. hard to say. Odd differences between grids. Many grids almost continuous deep wonky oscillations, also the grids with in general way higher amplitude than rest. If removed, not much is left, but could be salvagable? Sess3 the bad grids have huge line noise, separate ref maybe?
% R1159P (42) - 1sess - bad - fucked. Broad machine noise. Subsets of strips/grids having the same noise? Multiple amplifiers?? 
% R1080E (43) - 2sess - good-medium - sess1/2: finally again. Some stuff here and there but nothing major/continuous 
% -*R1142N (43) - 1sess - medium/good - some wonky oscillations, some spikelets, some oddly similar channels, but in general looks okay
% -**R1059J (44) - 2sess - good - sess1/2: a little bit of wonky oscillations, but not very wonky
% R1135E (47) - 4sess - medium/bad - sess1/2/3/4: many major events, but, ~1200 encoding trials, so maybe a decent number survive (however, very low accuracy..., ~10-20%). Also, few surface electrodes, and suspicious naming scheme suggesting hardware bipolars?
% R1147P (47) - 3sess - bad - sess1/2/3: fuckep up like R1159P, unfilterable noise all over the place but electrode set specifically, wth is going on. Haha in sess3 it even becomes phase-amplitude coupled....
% **R1020J (48) - 5sess (FR1/2) - good - sess1/2/3/4/5: some wonky oscillatory channels, but can be removed
% **R1045E (51) - 1sess - good - nothing major
%
% R1075J (50) - 2sess - 150:450 - allsurf - sess1/2: good, nothing much! USE
% R1084T (25) - 1sess -  53:247 - 
% R1120E (33) - 2sess - 207:393 - 40surf  - sess1/2: medium/good, prolly 30 surf remaining 
% R1129D (34) - 2sess -  40:188 - 
% R1151E (36) - 3sess - 208:548 - 20surf  - medium: maybe ~15 surf remaining BORDERLINE: check line noise and channels remaining
% R1155D (37) - 1sess -  38:87  - 
% R1166D (38) - 3sess - 129:771 - ~50surf - allsess - good (only 4Temporal)  VERY BORDERLINE: check line noise and channels remaining MANY LINE SPECTRA check whether 
%%%%%%%%%%%%%%%%%%%%%%% OLD NOTES %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visual artifact definition in samples (name is the ephys filebase, found in hdr.channelfile: hdr.channelfile{1}(1:end-4)
info.R1020J.artfctdef.visual.R1020J_26Jan15_1500       = [      1     20001     46674     91037     94053     97093    184948    189747    215198    249472    256779    264831    301444    306367    308001,...
                                                           433432    445093    499395    553872    572211    578210    749392    937118    948295   1023452   1037605   1048275   1153734   1178363   1218146,...
                                                          1235851   1237513   1334613   1428956   1541244   1542271   1739932   1758621   1770392   1840448   1938775   1940717   1941295   1944001   1944904,...
                                                          1948517   1949702   1951605   2040215   2189686   2260160   2282246   2301851   2310801   2323601   2336602   2340461   2349662   2478033   2506662,...
                                                          2525227   2553392   2556811   2573884   2598214   2633505   2650367   2653783   2794383   2822420   2936880   2949899   2963073   3068001   3166851,...
                                                          3229743   3242922   3349493;...   
                                                             1715     21840     47324     91920     95292     99058    186941    191622    219812    250606    258459    265804    302727    307211    310247,...    
                                                           434042    448403    500599    555013    573320    580316    750852    938945    950727   1024865   1038469   1050260   1156312   1179110   1219650,...
                                                          1236308   1238413   1337663   1430114   1542202   1543691   1740224   1759691   1773449   1843167   1940000   1941183   1943529   1944849   1947392,...
                                                          1949336   1950126   1953993   2042469   2190509   2260360   2283316   2303271   2311128   2326380   2340000   2343058   2350409   2480000   2507259,...
                                                          2526546   2554679   2557586   2575050   2599090   2633937   2652000   2655670   2794590   2823279   2938896   2950102   2966763   3068808   3167799,...
                                                          3232639   3243122   3350380]'; % mostly buzz, the occasional big event, quite some channels removed that were spikey/wonky (quite liberal). The occasional hint of mechanical artifacts
info.R1032D.artfctdef.visual.R1032D_FR1_0_08Mar15_1404 = [  13443     76801    105451    212854    296032    401928    411800    498483    815615    860582    890221    926124    959839    976548    992449,...
                                                          1034045   1104002   1230264   1234638   1259868   1261686   1283577   1367038   1395201   1611623   1634922   1787666   1984001   2004062   2006101,...
                                                          2071450   2079839   2102651   2140823   2218204   2909425   3061440   3250874   3282083   3454290   3531209   3564801   3723674   3813627   3914423,...
                                                          3919839   3930604   4040256   4074427   4082311;...   
                                                            15104     79229    106798    214097    297134    402389    412583    499200    816067    860828    891143    928000    960162    976893    992975,...    
                                                          1034867   1105198   1230584   1235200   1260239   1267899   1284596   1368157   1396128   1611865   1636020   1787990   1985646   2004330   2008815,...
                                                          2071787   2080162   2102846   2141843   2218501   2910239   3061760   3251848   3282290   3454957   3531533   3566149   3724290   3813873   3914914,...
                                                          3920162   3931611   4040582   4074841   4082712]';
info.R1033D.artfctdef.visual.xxx = []; 
info.R1034D.artfctdef.visual.R1034D_FR1_0_14Jun15_1542 = [  79808    239821    399666    443056    559593    576910    719673    879632   1039817   1199808   1359714   1519830   1544359   1556475   1679864   1839812,...
                                                          1999658   2159800   2319632   2323981   2462255   2479688   2557559   2639017   2679442   2799851;...
                                                            80304    240166    400364    444364    560282    579448    720252    880446   1040149   1200145   1360269   1520123   1546003   1560880   1680089   1840200,...
                                                          2000351   2160127   2320377   2325340   2464945   2480295   2559955   2645573   2681585   2800149]';
info.R1034D.artfctdef.visual.R1034D_FR1_1_15Jun15_1053 = [ 231920    239817    400943    519525    559628    719765    879834   1039830   1199795   1359778   1435808   1447365   1590544   1679821   1764337   1956943,...
                                                          1999709   2028801   2037507   2158256   2319825   2432001   2479834   2548630   2639898   2799752   2959692   3119821   3256199   3279671   3342363   3439782,...
                                                          3599731   3759739   3985498   4239722   4262401   4399757;...
                                                           232300    240123    402781    522540    560493    720304    880183   1040136   1200188   1360278   1438880   1451693   1594678   1680248   1767169   1958400,...    
                                                          2000274   2030700   2039280   2161357   2320149   2435293   2480157   2551190   2640175   2800274   2960269   3120265   3258915   3280222   3343134   3440136,...
                                                          3600489   3760295   3987200   4240248   4263844   4400136]'; 
info.R1034D.artfctdef.visual.R1034D_FR1_2_16Jun15_1325 = [79830    239709    399851    719825   1039795   1199808   1359855;...
                                                          80235    240291    400200    720162   1040218   1200170   1360235]';
info.R1045E.artfctdef.visual.R1045E_FR1_0_04May15_1104 = [ 171163    209565    301820    379019    384498    426380    430524    489693    573741    602921    668778    762922    777677    790086    878162,...
                                                           889369    960954    968491    983455   1053052   1077943   1117498   1197453   1271803   1337532   1354645   1461755   1559349   1611082   1635383,...
                                                          1639972   1658341   1754457   1776991   1848561   1853020   1955285   2021977   2042314   2047822   2121902   2160205   2204250   2241007   2317402,...
                                                          2346099   2355482   2360621   2365633   2410685   2433657   2437561   2444046   2447417   2464046   2500743   2503492   2547319   2597027;...    
                                                           172802    211615    303696    380723    386365    427211    432454    492163    575424    604593    670263    764103    778432    790893    879120,...  
                                                           891108    961611    969604    984755   1054306   1078920   1119819   1199318   1273163   1338166   1355578   1463467   1561205   1613403   1637682,...
                                                          1641187   1659408   1756451   1778103   1849070   1854144   1959085   2023742   2044074   2048105   2125872   2161836   2204822   2242682   2318991,...
                                                          2348935   2358910   2365448   2369436   2412071   2434973   2439200   2445344   2447670   2465171   2502519   2505168   2547597   2916214]';
info.R1059J.artfctdef.visual.R1059J_28Jun15_1200 = [   1    10299    21202    26351    53001    78875    82315    88981   102323   103577   144061   149440   244799   331175   378706   418617   462706   560457   670754   692001   734791   792872   968860  1049936  1093654  1172727  1300001  1334480  1492001  1660763  1667460  1688969  1711170  1732956  1789803  1970545  2062891  2100001  2227436  2336001  2476001  2852001;...
                                                    1961    12000    21844    26679    54780    79412    83336    90247   102868   103884   145534   150135   246840   332618   379046   418993   465139   562336   672000   692472   735679   794642   969203  1050046  1094380  1174453  1300530  1336000  1492970  1662485  1667654  1689268  1713332  1734586  1790961  1972000  2063558  2101711  2228000  2339324  2476832  2938000]'; 
info.R1059J.artfctdef.visual.R1059J_29Jun15_1039 = [   1     3145     8017    16815    26238    35053    36376    38202    41827    43871    47460    60001   167347   241932   451851   593726   733416   878912   887657   942690  1000872  1150807  1229771  1331383  1332775  1417819  1424001  1468936  1656900  1765730  1960775  1968001  2090158  2110142  2162383  2180307  2224727  2264791  2273928  2360771  2366960  2391827  2412001  2470246  2504001  2521259  2573186  2584001  2711512  2712457  2773916  2781408  2785311;...
                                                    1744     3433     8228    17296    28000    35154    36728    38501    43356    45381    57792    68000   168409   247658   452071   594143   735433   882336   888276   947755  1001562  1172000  1230497  1331707  1332929  1418340  1429647  1473288  1664000  1768000  1961284  1971118  2094074  2112518  2162546  2181804  2226235  2266405  2276252  2361889  2368000  2396865  2413997  2471122  2505389  2523715  2573941  2585324  2712000  2712824  2775582  2783275  2790000]';
info.R1068J.artfctdef.visual.xxx = []; 
info.R1080E.artfctdef.visual.R1080E_FR1_0_25Sep15_1253 = [ 141029    160546    336736    412262    429732    568303    608670    879121    956954   1069289,...
                                                          1098167   1242757   1316505   1342657   1360071   1387479   1408925   1514485   1658341   1683984,...
                                                          1700286   1793471   1914528   1963326   2052293   2086424   2135601   2247698   2368404   2377621,...
                                                          2429157   2533674   2545453   2571120;...
                                                           141646    163648    337325    416963    432069    572542    611388    882090    957559   1070928,...
                                                          1100404   1245738   1316977   1344168   1361542   1389590   1411975   1517313   1661358   1686652,...
                                                          1704106   1795736   1917549   1965381   2056310   2088487   2140786   2251300   2371023   2381616,...
                                                          2431777   2536269   2550316   2572066]';
info.R1080E.artfctdef.visual.R1080E_FR1_1_01Oct15_1149 = [  20915    101830    180592    198274    309276    488351    504230    715285    833646    861050,...
                                                           943524   1039517   1153294   1272855   1372723   1374899   1537864   1646752   1691823   1847989,...
                                                          1927893   2047681   2070110;...
                                                            23976    104812    181958    199800    311688    490192    507492    717424    835164    863758,...
                                                           947690   1042792   1155426   1274724   1374624   1377183   1540274   1647707   1695063   1848329,...
                                                          1928982   2048153   2071076]';
info.R1128E.artfctdef.visual.R1128E_FR1_0_11Jan16_1359 = [  48233     76895    118006    164428    213524    348612    377140    646724    648607    671490    847943    956817    960611   1056226   1067825,...
                                                          1119373   1120511   1122084   1147769   1239915   1263981   1264412   1306693   1352032   1355039   1379110   1390386   1390838   1546980   1624557,...
                                                          1670106   1703003   1710009   1746430   1819228   1832686   1832978   1939874   2046557   2146507   2146586   2346706   2435069   2695429   2745978,...
                                                          2825866   2845953   2945921   2950190   2965785;...   
                                                            48838     77363    124327    165369    214475    350747    377632    647352    650243    672490    850194    958153    961762   1057314   1068808,...   
                                                          1119814   1121393   1122257   1149988   1239969   1264031   1264467   1307102   1353566   1355092   1379487   1390610   1391400   1549751   1625458,...
                                                          1670196   1703724   1710564   1746480   1819965   1832838   1833028   1940374   2049203   2146557   2149030   2348465   2436141   2695495   2748584,...
                                                          2827306   2848043   2948260   2950436   2966194]';
info.R1135E.artfctdef.visual.R1135E_FR1_0_29Jan16_1259 = [18482    31582    57671    67441   114428   135983   179821   192792   201901   207191   215298   242817   259741   278211   283967   300840   305775   312961   319931   322855   361628   410952   443971   451549   467845   471287   477579   498126   501582   504263   515730   554464   582509   587413   594878   603397   608875   611389   644149   649877   655823   661288   673963   688771   700705   707293   711743   757417   813869   852878   866832   886395   900323   909971   919081   990721  1014985  1022977  1031804  1056905  1065021  1068936  1078612  1083162  1090149  1098183  1105225  1112653  1118606  1123946  1165646  1170093  1178297  1207725  1230317  1268532  1270839  1283821  1291453  1294877  1303967  1306518  1311763  1320708  1324119  1370065  1374050  1381540  1459723  1466533  1478102  1480377  1493543  1501425  1506493  1518016  1523068  1548333  1574145  1578421  1649326  1665570  1693692  1741411  1750249  1759103  1858660  1885712  1898101  1968468  2036172  2061571  2082223  2095787  2137149  2162237  2168260  2248009  2252442  2263531  2273725  2280720  2287313  2296618  2354402  2372832  2377196  2400398  2600285  2672932  2679077  2701297  2707669  2712347  2717281  2720492  2749768  2760342  2764738  2768455  2773225  2786669  2792759  2797201  2801197  2804470  2814316  2823229  2825455  2842277  2859673  2873402  2876323  2895736  2899891  2907636  2917081  2930567  2937061  2945977  2957229  2959758  2963835  2970160  2973025  2981774  2991281;...
                                                          19980    37319    59440    72653   115884   136696   181351   193899   203796   211788   216047   244540   260956   279720   285294   302237   307363   314370   321025   323676   364665   416652   446136   452995   469164   472202   478485   499423   503496   505743   516788   556273   584189   589315   596698   604510   610157   616321   647352   651348   657014   663336   675324   690534   702097   709891   714158   758304   815184   853827   868759   891108   901645   911088   920427   992299  1019380  1024831  1032992  1058208  1066021  1070928  1080745  1084497  1091549  1099719  1106407  1114571  1120176  1125928  1166832  1170828  1179875  1209597  1231879  1269782  1272266  1285908  1292673  1296650  1305936  1308139  1314684  1321810  1325803  1371412  1375225  1383271  1460916  1467560  1478780  1481664  1498500  1503129  1508484  1519282  1524578  1548548  1575594  1580398  1650688  1666332  1694873  1742256  1751297  1760900  1859557  1887776  1898940  1970028  2037062  2062669  2083433  2097244  2137860  2163551  2169308  2249392  2253744  2270944  2276135  2281716  2289038  2297700  2356087  2373624  2378769  2401890  2601336  2673742  2680619  2702633  2708720  2714001  2718657  2722521  2751216  2761695  2765914  2769228  2774902  2789097  2793853  2798864  2802836  2805941  2815783  2824806  2826603  2843497  2861049  2874871  2877120  2898023  2901982  2910336  2918371  2931957  2938577  2947022  2958769  2961036  2965824  2972416  2979057  2982951  2992136]'; 
info.R1135E.artfctdef.visual.R1135E_FR1_1_30Jan16_1457 = [      1      6843     10555     26568     49049     55692     63937     67503     71929     76575    112547    117075    119188    125051    129484    132208,...
                                                           146236    153447    160988    167258    184067    188106    191610    208591    213902    222034    232363    237994    240293    247753    257286    260246,...
                                                           278179    286558    334464    342523    346893    354823    361709    367633    386450    393328    400482    440136    443896    446046    456383    459541,...
                                                           468516    480158    486218    489632    494396    502538    507493    530113    538188    541883    553146    556683    575425    590560    600456    612882,...
                                                           640430    651349    664027    668821    671329    685854    688130    700190    702601    712430    753910    761325    764564    775725    779396    781954,...    
                                                           787713    799201    803874    811887    815185    819181    865540    868398    871401    880346    892763    903417    907093    914107    917657    921688,...    
                                                           925580    969928    987408    989271   1013086   1015689   1021207   1027521   1048876   1085965   1118507   1121569   1132652   1146853   1186813   1191169,...
                                                          1207489   1210789   1236720   1238761   1243233   1256517   1261960   1298346   1305677   1309341   1313498   1384454   1399079   1437081   1444890   1451803,...    
                                                          1455896   1461760   1478521   1484585   1490320   1498687   1544178   1546494   1553521   1555683   1564502   1586413   1592337   1599690   1660127   1664241,...    
                                                          1672437   1675155   1689165   1696064   1699724   1707249   1713823   1719839   1758910   1767691   1769820   1787072   1789647   1796708   1799477   1802850,...
                                                          1815783   1824940   1858627   1873499   1877766   1886540   1897193   1899755   1915275   1922335   1924550   1963651   1984920   1992310   1994462   2001997,...   
                                                          2015859   2044725   2066728   2083867   2087669   2093905   2113885   2126805   2146304   2183388   2191176   2215528   2235637   2280726   2284706   2292222,...   
                                                          2301138   2305526   2322579   2335776   2338400   2355611   2380032   2412116   2456523   2506216   2522044   2553227   2568000   2612348   2621896   2629103,...
                                                          2644013   2654615   2658759   2665333   2674133   2706893   2720723   2722633   2724260   2730937   2741706   2744170   2749249   2772027   2810806   2824442,...
                                                          2845153   2863986;...   
                                                             1347      8529     11988     27972     50913     60759     65206     68965     73955     78034    113952    117750    121099    126703    131346    135310,...
                                                           148171    154708    163279    169142    185945    189857    193080    209907    215015    223699    236086    239760    241236    250112    258329    261533,...
                                                           279999    296152    336528    344649    348930    356361    364063    368767    387829    395929    401764    441270    444847    450400    457678    461894,...
                                                           470544    483953    487512    490618    496188    504178    508394    531456    539822    542834    554412    557946    576659    591408    602168    614000,...
                                                           642026    653802    665330    670196    673043    687312    697999    701727    710232    713930    755108    763236    765510    776904    780218    783876,...
                                                           789733    799661    804774    813102    816432    820672    866905    869478    872508    882137    894359    904790    910265    916131    919080    923076,...
                                                           927072    971608    988349    991835   1014422   1019847   1022577   1028886   1049016   1087760   1119277   1123439   1133987   1148089   1188004   1192902,...
                                                          1208808   1211641   1237508   1240199   1245626   1258740   1263189   1299457   1306928   1310688   1314605   1385813   1400299   1438051   1445683   1452900,...
                                                          1457218   1462536   1480047   1487006   1491480   1500183   1545004   1548524   1554444   1556393   1565835   1588498   1593758   1601023   1661489   1665375,...
                                                          1673208   1679012   1690308   1698300   1700775   1708544   1714821   1720860   1760036   1768804   1770934   1787887   1790755   1798118   1801478   1804542,...
                                                          1817494   1826117   1859796   1875117   1878879   1888423   1898693   1901150   1916943   1924049   1926072   1964981   1986793   1993455   1995902   2002704,...
                                                          2017181   2045945   2068042   2084348   2088443   2095517   2115092   2127926   2147777   2185166   2192694   2216624   2237184   2281621   2285712   2293378,...
                                                          2302794   2306435   2323617   2336883   2339762   2356624   2381440   2413584   2457540   2507309   2524101   2554663   2569394   2614181   2622946   2630836,...
                                                          2645352   2656326   2660288   2666395   2675466   2709174   2721694   2723579   2727583   2732159   2742958   2745252   2750367   2773742   2811486   2826216,...
                                                          2846679   2865044]';
info.R1135E.artfctdef.visual.R1135E_FR1_2_01Feb16_1500 = [75925    79921   238055   243757   255745   262364   266570   271022   276780   311689   348695   353572   362962   419320   427678   444940   450244   457862   511220   530411   573129   577718   616698   651878   657727   697265   706667   708590   711289   734854   753848   758878   774760   777666   831489   875619   925147   942809   961659   986701  1063332  1068375  1081845  1090575  1127499  1138471  1149670  1216595  1248630  1275049  1424677  1493028  1505030  1525820  1548164  1550293  1555151  1558441  1597241  1601953  1615347  1641468  1651643  1678321  1709730  1739249  1750249  1772420  1777930  1780383  1787067  1793912  1848371  1855976  1893847  1943432  1954045  2064359  2070517  2092790  2105893  2112942  2115880  2145853  2155513  2163800  2186162  2210726  2215632  2245753  2255882  2260517  2285460  2301697  2309498  2314676  2365415  2369408  2404492  2406761  2409288  2416971  2423956  2450980  2455282  2460103  2467507  2476441  2477548  2507639  2529469  2566795  2583310  2603014  2608000  2611322  2618549  2622365  2653218  2679665  2688718  2699661  2703869  2717047  2722343  2725939  2762266  2771170  2809302  2812229  2815129  2863516  2920099  2924130  2963556  2989785  3008570  3022908  3060937  3064933  3088191  3105079;...
                                                          77416    84273   238916   244929   256724   263736   267732   272498   278245   318557   349896   355098   365064   423576   429478   447552   451548   459190   511770   539460   575007   579420   618866   654094   660709   698504   707292   710371   711983   736077   756304   759484   775954   779220   833172   876643   927072   945504   962543   988196  1064608  1069482  1087255  1091469  1129242  1139601  1150726  1217463  1248876  1276393  1425486  1497218  1507024  1527068  1549056  1550858  1556130  1559361  1598334  1603247  1616631  1642595  1656081  1679090  1714284  1747642  1755194  1773471  1779310  1781573  1789431  1805015  1850148  1867882  1896125  1945539  1955088  2066310  2072978  2094336  2107301  2113884  2117880  2147506  2157840  2167134  2187189  2213034  2221203  2247623  2257740  2261736  2288367  2304276  2310837  2317680  2366287  2370540  2405592  2407831  2410358  2418525  2425572  2452313  2457149  2461536  2470258  2477520  2480218  2509081  2530445  2568541  2584141  2604180  2610182  2612293  2619973  2624388  2654192  2680697  2689740  2701296  2708881  2717562  2724390  2727108  2765232  2773224  2810670  2813184  2817180  2865132  2921457  2925469  2965032  2991000  3014554  3024154  3062302  3068647  3090137  3106138]'; 
info.R1135E.artfctdef.visual.R1135E_FR1_3_02Feb16_1126 = [35145    39961   103897   108935   117864   162121   219781   238500   240932   247767   264578   343119   360449   367633   428282   430760   435925   446333   457910   481304   489774   517137   537804   559960   563162   571429   586129   588641   622992   635365   661452   715707   735033   751865   767233   828411   835165   843157   868774   919208   947053   948868   967033  1012270  1037441  1043808  1050949  1062246  1073123  1086451  1100961  1154845  1166833  1176112  1249709  1259347  1266733  1287785  1295513  1351517  1373843  1402597  1419540  1448057  1459084  1470529  1494991  1511362  1519247  1553932  1563184  1591021  1604642  1611861  1645166  1666333  1717252  1753616  1761863  1845709  1849483  1858141  1905115  1908142  1948164  1950484  1973350  1984130  1989265  2005329  2046767  2053945  2088139  2091816  2109351  2141601  2153240  2162659  2220074  2263415  2269729  2301697  2305693  2321677  2326629  2362048  2370446  2373625  2378483  2413585  2421910  2468658  2481517  2495906  2499630  2518005  2565290  2577071  2583780  2612442  2659188  2675312  2693883;...
                                                          37075    40985   105023   111041   118926   162911   221377   239596   242624   251748   266187   351432   361587   369449   429698   432333   438181   448118   458752   486987   491021   529802   542821   561405   566794   573055   587412   589773   623376   639360   662720   717550   738620   753664   767991   830138   838195   849812   869677   921335   948064   949208   967942  1014129  1041011  1047969  1056020  1063142  1074337  1087524  1104520  1162202  1170080  1177413  1254063  1263858  1267798  1290708  1298191  1353445  1375679  1405670  1422576  1456705  1460363  1477668  1496421  1514484  1521495  1556057  1565461  1593232  1609114  1613282  1647571  1668998  1719184  1755825  1765208  1846756  1851176  1874124  1906092  1909378  1949382  1951975  1974287  1986012  1991170  2007641  2049718  2055640  2089324  2092998  2111842  2143456  2155089  2164392  2221945  2264576  2270630  2303038  2310423  2323246  2329315  2363198  2372251  2377620  2379469  2414808  2423396  2469912  2483212  2497902  2500901  2525472  2567290  2579192  2585209  2614405  2660105  2677811  2695112]';
info.R1142N.artfctdef.visual.R1142N_FR1_0_22Feb16_1029 = [   1     6609    14944    35637    45472    46716    49932    64344    78538    81908    88001   100001   148001   156001   311063   312001   332001   348001   364001   372001   382971   400001   412428   420001   428880   433551   481263   488050   495891   504301   512001   517307   645335   652001   712001   737130   822293   851174   900973   940001   972783  1032803  1160632  1254263  1478045  1629835  2192001  2218754  2557327  2749831  2818145  2849795  2873609  2926944  2946307  3233577  3239049  3254686  3272001  3300767  3355669  3621311  3675077  3764001  3793013  3965247  4020912  4034779;...
                                                          3175     9211    16000    36000    46389    46874    50554    64607    81798    85248    93034   144000   152000   300000   311600   324000   336000   352000   368000   372744   388000   404000   412818   428000   430425   433654   482917   492000   496010   505966   516000   519167   647973   656000   716000   738542   832000   852000   901820   945429   973369  1033473  1160826  1255489  1478332  1631013  2192647  2219025  2564595  2750139  2819653  2850276  2874368  2931864  2949127  3233897  3239969  3261219  3272546  3301727  3356599  3622264  3681320  3764812  3794094  3965973  4021711  4035650]';
info.R1147P.artfctdef.visual.R1147P_FR1_0_12Mar16_1734 = [    157     13332     18950     24001     30407     42511     48394     61990     80001     99479    137178    154514    180001    208585    234622,...
                                                           246597    335715    357888    365880    379138    434638    470611    479765    498614    520692    539022    548001    570963    578768    668001,...
                                                           707248    765810    784001    823449    858955    876278    892001    895374    904227    984001    996001   1028213   1079315   1108458   1112001,...
                                                          1126399   1215218   1336552   1340585   1437825   1489095   1544001   1621157   1652417   1757423   1768388   1776001   1885187   1936305   1984434,...
                                                          1988001   2098495   2141708   2205426   2310012   2410025   2420001   2428011   2461800   2489165   2525601   2529106   2611701   2621485   2681694,...   
                                                          2716001   2739380   2763734   2924958   2953154   2967673   2987103;...   
                                                             1360     14198     22526     25510     31754     48000     50281     63633     82013    101029    138865    156000    181679    211160    238214,...    
                                                           252000    337360    358816    368927    380739    436212    471477    481727    500000    521406    539867    549007    572000    579460    671190,...
                                                           708586    766910    785365    824000    860000    877513    893179    897306    905349    986905    998590   1028462   1080448   1109618   1116000,...
                                                          1126985   1218641   1337013   1344000   1440816   1490811   1552000   1621819   1654787   1760000   1768895   1778440   1887999   1937411   1987953,...
                                                          1990867   2101620   2142435   2208000   2311284   2411574   2420491   2429423   2462925   2489548   2528743   2531448   2615800   2622897   2682150,...
                                                          2720000   2740647   2764956   2926792   2953870   2969044   2988000]';
info.R1147P.artfctdef.visual.R1147P_FR1_1_20Mar16_1442 = [   2842     18880     22532     26646     33759     38460     44001     53340     56706     57875     60069     63879     65751     70638     82686,...
                                                            93041    107113    159001    162950    172001    176845    181152    184001    189568    193684    202138    205264    212001    289928    292560,...
                                                           330947    335576    343444    376184    390815    396557    494729    592824    608001    630351    648001    697439    744001    795130    805493,...
                                                           838436    841651    871066    947105    950702    960001   1052746   1060001   1161450   1264834   1269192   1280001   1304001   1330630   1354189,...
                                                          1357009   1395250   1413625   1422753;...    
                                                             4000     21518     24600     29416     35934     41607     51391     55714     56944     59488     62160     65229     67566     71491     85669,...
                                                            95413    109303    160000    169427    174991    179808    182660    185577    191343    195722    204054    208000    216000    292000    294789,...
                                                           332000    336827    344620    376582    394332    397652    497037    596000    609943    631400    650300    699195    745317    795948    807020,...
                                                           839561    842789    871654    949007    952000    961005   1054800   1062201   1168663   1268000   1272000   1282306   1307964   1332000   1356586,...
                                                          1364000   1400000   1416000   1424000]';
info.R1147P.artfctdef.visual.R1147P_FR1_2_22Mar16_1329 = [     1      8969     12001     21450     32001    101380    112770    133885    151221    181591    228757    236558    246044    252001    260001,...
                                                          274694    360267    376001    432001    437740    452001    496001    503000    544001    564746    588963    669377    768544    799506    856001,...
                                                          881748    889891    892869    901842    918377    970090   1018235   1085563   1089418   1109995   1296705   1300910   1309471   1332001   1391215,...
                                                         1529700   1532488   1547853   1636001   1714140   1734103   1772364   1783220;...    
                                                            5378     11149     14273     23077     33650    104000    118644    135421    152297    182019    231590    239980    247222    257362    261198,...
                                                          276000    363375    377857    433693    439063    455295    499031    503885    545470    565988    590639    671276    770074    809862    857953,...
                                                          883880    891082    894566    913674    919596    971359   1019880   1087663   1091773   1112000   1298706   1302521   1310604   1334273   1392365,...
                                                         1531445   1533292   1548363   1638297   1715282   1735112   1773994   1786582]';
info.R1149N.artfctdef.visual.R1149N_FR1_0_07Mar16_1014 = [   1     4001     9521    22557    40001    49819    60001    80001    90464   134859   138581   143278   173702   193408   208888   214198   216001   224001   232001   235710   244001   266847   270666   274766   277085   288001   291488   316900   372001   384001   410412   424001   468001   496658   504001   512001   528001   596001   616001   625984   636001   648001   660001   669037   676001   682057   696001   700001   708001   717521   736001   808001   857114   900001   941033  1054904  1064936  1090142  1113251  1123234  1146053  1177553  1185110  1200001  1220001  1225251  1247460  1274166  1278392  1316001  1325065  1346053  1380896  1426162  1554105  1576001  1649573  1652001  1656001  1673867  1678468  1684001  1693182  1752985  1766117  1828521  1850202  1948001  1960001  1968283  1977734  1995153  2021017  2042178  2044001  2049755  2099710  2126081  2136446  2237674  2240094  2249888  2278795  2281936  2288001  2296001  2348001  2356001  2400001  2408295  2495121  2524860  2549694  2555182  2568497  2676565  2690021  2694754  2705327  2709674  2716332  2729190  2736840  2740908  2747214  2761384  2768001  2772001  2784001  2811375  2819041  2828001  2836432  2844190  2853227  2860719  2875286  2883367  2896884  2902988  2910246  2918142  2931363  2943219  2949694  2991738  2998488  3005037  3011661  3017134  3028094  3035407  3049460  3058803  3061448  3084001  3091669  3099315  3112960  3123012  3133839  3149146  3158263  3160985  3190202  3219291  3232598  3253235  3262811  3268001  3276332  3283037  3285287  3296880  3316182  3330367  3349331  3362480  3377384  3388001  3398343  3406662  3420001  3456642  3464001  3487250  3507264  3519403  3524420  3533747  3549952  3557142  3596001  3606436  3615101;...
                                                          3433     6522    18364    23340    48000    50167    60963    84000    98223   135525   139025   144000   175352   195416   209812   214945   220000   225240   234199   236679   246905   268000   272000   276000   279098   288643   292369   318594   372474   384615   411537   431328   469449   498280   505183   524724   532000   600000   620000   628000   645369   656000   662433   670574   678050   684000   697171   704000   714961   728000   740000   809022   863473   902481   942909  1057252  1067513  1092635  1114457  1124000  1150042  1178969  1186594  1201877  1224000  1227106  1251860  1276000  1279771  1320000  1326114  1352000  1384000  1427207  1556000  1584000  1650348  1653735  1660000  1676000  1680377  1686981  1694143  1755481  1767070  1830126  1851674  1953812  1963860  1970332  1980000  1996679  2024000  2043860  2045445  2050328  2102268  2127904  2140953  2238649  2242167  2250187  2280000  2284000  2292000  2300000  2350675  2359062  2404000  2412000  2499259  2527590  2552000  2556845  2570707  2680000  2692000  2697147  2707602  2712000  2722284  2730977  2738348  2744000  2751009  2764000  2770002  2773453  2786497  2813143  2820990  2829683  2838518  2850251  2855090  2862743  2880000  2888349  2898445  2904385  2913881  2920000  2932599  2944031  2954780  2992982  3000000  3006239  3012716  3018832  3029832  3044000  3052000  3059844  3063013  3087142  3093469  3101481  3118739  3124000  3135058  3150292  3159505  3163404  3192578  3220849  3234082  3255102  3264603  3269788  3279731  3284000  3289050  3299634  3322501  3331283  3351058  3363771  3380000  3388748  3400000  3412000  3424000  3458594  3466570  3489905  3509449  3520897  3529953  3534743  3550147  3559348  3596615  3607997  3616000]';
info.R1154D.artfctdef.visual.R1154D_FR1_0_10May16_1745 = [20235    37626    70279   157493   162968   175675   218081   234134   326809   428001   490178   618476   636001   924001  1020001  1040154  1048001  1104739  1317783  1321017  1417904  1537924  1740682  1826166  1845589  2126367  2331194  2535198  2639089  2774025  2798581;...
                                                          23501    37886    72000   159017   163586   175786   220000   235537   328000   432000   496000   621324   640000   926671  1021703  1040986  1052000  1107554  1319574  1324000  1419449  1540000  1741123  1829465  1847445  2133498  2334530  2541977  2642872  2775493  2808000]';
info.R1154D.artfctdef.visual.R1154D_FR1_1_16May16_1545 = [20001    27123    67381    74353   104954   116594   124245   140001   164852   186263   208116   225618   240870   258621   279313   334589   348829   372833   384001   429753   533218   540549   550140   587115   642218   649412   675490   691349   748001   800561   875452   883819   927542   956001  1001351  1029136  1035411  1058589  1065039  1084706  1158775  1165940  1176731  1183266  1193456  1253984  1260001  1266226  1281170  1357513  1364134  1410327  1412569  1466525  1565880  1586908  1629138  1663581  1672767  1744001  1750396  1767524  1789307  1872247  1881493  1889093  1969279  2053722  2090682  2116767  2160001  2222799  2266565  2272763  2317013  2381392  2464670  2573154  2596755  2675762;...
                                                          22622    32456    70451    77030   106624   117997   125385   141097   166310   187733   209080   231326   243346   260000   281314   337915   350177   374124   388000   436490   535275   544768   551749   588929   646911   651372   676818   693919   753300   803239   876843   886497   928526   959158  1003364  1030511  1037217  1062290  1067499  1088506  1162292  1171013  1180000  1184365  1193887  1256000  1262848  1268000  1282913  1358651  1364945  1411445  1414933  1468000  1570808  1588845  1630401  1665481  1674179  1746518  1752000  1782179  1790441  1874719  1883336  1892000  1972000  2056728  2093373  2118211  2161985  2225260  2269852  2279364  2318767  2400000  2478187  2574489  2599384  2908037]';
info.R1154D.artfctdef.visual.R1154D_FR1_2_17May16_1818 = [     1   313904   334702   345041   347000   364767   432384   450396   530924   548001   554899   564283   633493   641763   660916   739774   757726   768892   846875   850468   857509   863262   943734   974029  1039295  1073956  1098670  1254097  1342988  1359532  1369146  1452001  1530545  1533476  1544557  1552360  1556291  1573396  1578589  1643097  1700001  1761533  1826283  1843230  1856001  1868787  1944001  1956960  1970496  1994851  2040001  2053315  2075815  2082881  2161162  2168360  2178021  2184819  2196102  2262819  2336698  2375121  2410416  2419714  2456662  2487383  2501452  2508001  2514779  2549166  2562448  2580606  2594069  2604610  2613448  2644001  2649880  2674009  2711157  2787065  2802525  2893448  2904190  2973299  2976001  2987778  3034383  3096783  3109791  3133227  3191093  3195682;...
                                                          251755   315727   336740   345340   347392   367287   435283   452389   544000   551336   556000   566264   640000   644861   664820   744655   758312   770106   848663   852929   858640   866122   948000   975094  1043838  1076139  1100000  1255574  1348635  1362788  1370747  1453498  1532000  1533893  1547868  1553691  1557788  1573554  1579582  1644361  1701518  1762792  1828788  1845187  1857885  1871537  1951110  1960000  1972329  1996000  2043090  2055658  2076587  2084458  2162800  2172000  2179086  2185695  2197848  2266239  2338981  2378227  2412000  2421099  2458453  2490082  2501977  2509530  2516000  2551642  2564744  2581808  2595080  2606380  2614046  2644607  2651066  2676000  2712000  2788000  2804000  2895501  2904744  2975844  2979759  2989481  3035932  3099150  3110570  3134425  3191541  3196171]';
info.R1156D.artfctdef.visual.xxx = []; 
info.R1159P.artfctdef.visual.R1159P_FR1_0_02Apr16_2014 = [121227    576594    650992    714565    719210    830134   1144920   1386847   1578166;...
                                                          121760    576893    653107    715421    719465    830731   1145183   1387247   157874]'; 
info.R1162N.artfctdef.visual.R1162N_FR1_0_18Apr16_1033 = [   1    30723    37480    50108    61501    69919    84481    91438   101773   119699   131747   164696   175692   183908   190910   199730   202215   211267   219387   227948   237801   261472   270778   275690   282068   288001   296799   322664   328671   344892   356337   359360   363441   419156   438889   448394   463301   532367   700267   764001   777958   929020  1062055  1132001  1139245  1151789  1285582  1425872  1522114  1542590  1618090  1660980  1724189  2023019  2105993  2129348  2231084  2416539  2486796  2584001  2595272  2614772  2692001  2708001  2712001  2717923  2725751  2749810  2973485  3000001  3023001;...
                                                          1521    31317    38269    50673    64960    72553    85019    92098   102500   122002   132448   165369   176809   184500   192139   200111   202902   211794   219999   228497   238667   261889   272641   276308   282764   290472   298009   323348   329240   345490   357052   360030   364852   420602   439204   448395   464074   534354   702206   765720   787192   930172  1064139  1136000  1141373  1152874  1287090  1428000  1524231  1543558  1620000  1661760  1727628  2024000  2107329  2131136  2235012  2420000  2488513  2588000  2596000  2623088  2694023  2710276  2713408  2719571  2727039  2750236  2974569  3004792  3025730]'; 
info.R1167M.artfctdef.visual.R1167M_FR1_0_25Apr16_1912 = [ 60259    169920    280912    319855    388783    407085    522480    620001    643561    727512,... 
                                                          828001    907774    959839   1083718   1223420   1496275   1707375   1826416;...
                                                           63267    172651    283892    320651    390268    408000    525103    622070    645816    728433,...
                                                          828691    908651    961582   1084349   1224308   1497691   1707719   1828502]';
info.R1167M.artfctdef.visual.R1167M_FR1_1_28Apr16_0038 = [165686    393146    409734    433996    463012    497537    509102    543641    589968    777271,...
                                                          849255    938029   1139004   1348001   1400001   1654396   1730259   1821533     12646;...
                                                          166755    394086    410586    435142    463687    497856    510824    544000    591408    777707,...
                                                          851441    938429   1140639   1348639   1400587   1655078   1730634   1823110     14722]';
info.R1175N.artfctdef.visual.R1175N_FR1_0_19May16_1858 = [      1     30988     65601     79995     85558    126610    223040    248167    311973    315701    323179    324210    328001    332101    345368,...
                                                           350884    353012    357974    360801    387781    389639    468605    522062    537101    541201    608195    635501    638540    647310    708344,...
                                                           716426    760680    818240    824776    877401    889701    895080    897901    914759    926965    996301   1024648   1068825   1211598   1264396,...
                                                          1294779   1313285   1446378   1451401   1463701   1470256   1516684   1531761   1547935   1552906   1553901   1568573   1628891   1637108   1652301,...
                                                          1678538   1682439   1726884   1733311   1749937   1756980   1763572   1767101   1785791   1802323   1808101   1813576   1818208   1825460   1833624,...   
                                                          1849614   1857301   1861401   1866146   1881901   1948080   1972619   2057457   2066401   2073427   2089576   2094668   2106905   2167245   2175417,...
                                                          2182559   2187742   2190385   2193501   2199191   2205481   2209901   2283199   2287062   2291901   2312401   2404061   2444279   2490269   2513883,...   
                                                          2518732   2528783   2559087   2566601   2590509   2621678   2625938   2644567   2655993   2656801   2665001   2672707   2673201   2686986   2689601,...   
                                                          2706698   2721304   2727165   2730783   2736434   2741473   2742901   2767501   2781878   2786372   2798215   2806746   2814373   2816701   2876815,...
                                                          2880444   2882301   2886401   2894601   2904344   2908794   2911001   2920036   2933661   2942990   2959994   2987135   2989995   3011503   3054501,...   
                                                          3070380   3076442   3092269   3097774   3103001   3108303   3194984   3206201   3210301   3218501   3235714   3239001   3243672   3258991   3262694,...   
                                                          3264555   3266725   3328258   3333301   3376034   3401317   3441967   3444001   3453430   3483929   3556167   3564871   3601027   3634893   3646598,...
                                                          3654005   3657201   3669939   3673601   3741686   3798220   3838701   3892697   3897370   3996283   4028110   4036627   4044929   4092317   4109568,...
                                                          4120501   4125020   4186101   4201170   4203080   4220987   4239401   4249122   4264112   4282586   4303926   4306555   4313201   4367375   4409234,...   
                                                          4452601   4460801   4469001   4473101   4493601;...   
                                                             3848     32388     67834     81613     86100    127100    224416    249019    315354    317912    323900    327817    332013    340037    347730,...    
                                                           351447    355941    359036    364879    388413    391544    468926    523498    539623    543241    608733    636257    639287    647642    709592,...
                                                           717289    761221    818643    825160    885139    890689    895422    902000    915063    927401   1003240   1025377   1069093   1212651   1265265,...
                                                          1295411   1313763   1449316   1455080   1467800   1470766   1517000   1532994   1549044   1553717   1556861   1578500   1629600   1637399   1662601,...
                                                          1678886   1686615   1727356   1734300   1750700   1757547   1764242   1771200   1785995   1806523   1809601   1814007   1824500   1825692   1835501,...
                                                          1849868   1858980   1865500   1866389   1882258   1948635   1981432   2058756   2067433   2084450   2089847   2094831   2107134   2167485   2175619,...
                                                          2183048   2187936   2190612   2197600   2199999   2205799   2217104   2283700   2287800   2300100   2312708   2407355   2444495   2490540   2516184,...
                                                          2518984   2528991   2559581   2574800   2599400   2621957   2626035   2650748   2656336   2663178   2667413   2673050   2685500   2687907   2696261,...
                                                          2714200   2721762   2728913   2734700   2738250   2742791   2763400   2778003   2782152   2792524   2798525   2810318   2815875   2818154   2877656,...
                                                          2880712   2884545   2890750   2895861   2906900   2909205   2916276   2921825   2935600   2943680   2961921   2988900   2990247   3012829   3055337,...
                                                          3070626   3080254   3092995   3099600   3103700   3108505   3195569   3207456   3215187   3219456   3236101   3241341   3252390   3259435   3263042,...
                                                          3264776   3267013   3331621   3337400   3376390   3401792   3442423   3447413   3476800   3484114   3558800   3565238   3610301   3636700   3651418,...
                                                          3655116   3661300   3670535   3726900   3742938   3799343   3838853   3892919   3897633   3996463   4028541   4037749   4045404   4097095   4111293,...
                                                          4122816   4125277   4191505   4201468   4203395   4235300   4240579   4249765   4264785   4283080   4304544   4306834   4316337   4367742   4411600,...
                                                          4455347   4464900   4470889   4489500   4571000]';
%
info.R1075J.artfctdef.visual.R1075J_22Aug15_1233       = [      1     59081    110174    160017    219065    249399    259877    348001    349845    372001    474573    499669    506597    660336    849839,...
                                                          1337541   1591894   1774866   1776001   1944165   2336001   2345726   2361237   2382416   2412001   2686912;...
                                                             1546     61135    110348    163521    220260    251042    260262    349038    351222    372744    475902    501749    506797    660755    850159,...
                                                          1338380   1592094   1775314   1776804   1946932   2337467   2348000   2361437   2384986   2412986   2688000]';
info.R1075J.artfctdef.visual.R1075J_25Aug15_1409       = [       1    167095    200001    264803    277439    281722    313485    366085    371608    398762    572001    591234    697025    754069    756207,...
                                                            848557    856186    947811   1149896   1198525   1270339   1390029   1460001   1490605   1846666   1968348   2042531   2060001   2163054   2168150,...
                                                           2212001   2252001   2304001   2498379   2524231   2568001   2606359   2728811   2741436;...
                                                              1425    168000    203571    266497    278451    281989    314771    366711    373050    399388    573198    592000    697582    754485    756933,...   
                                                            849373    857433    948413   1150139   1201361   1271227   1396300   1461272   1496000   1848000   1976000   2042581   2068000   2164000   2175205,...
                                                           2215570   2258832   2305707   2500000   2525381   2572000   2610542   2730739   2745981]';
info.R1120E.artfctdef.visual.R1120E_FR1_0_18Dec15_1403 = [19167     68912     85665    192932    657858    679321    680817    979623   1328066   1362826   1403753   1998480   2313983   2697889   2714469;...
                                                          19980     70278     86367    194184    658797    680549    680818    980899   1328929   1364486   1404938   1999182   2318991   2699038   2715104]';
info.R1120E.artfctdef.visual.R1120E_FR1_1_19Dec15_1047 = [  49648     53483    103897    149488    332728    460250    551807    749142   1102897   1148637   1310689   1421416   1481139   1610389   2013585,...
                                                          2017981   2317681   2325673   2347493   2417581   2641736   2665583;...  
                                                            50121     57058    104941    150109    334678    461097    553419    749364   1113314   1149158   1318680   1422130   1482223   1614384   2015146,...  
                                                          2025972   2321676   2329668   2347881   2419511   2645352   2667028]';
info.R1151E.artfctdef.visual.R1151E_FR1_0_11Mar16_1849 = [   1     18393     34318     38662     52579     55183     74229    135695    150397    160001    187116    191511    211573    256262    278633,...
                                                        319346    321490    332001    368001    376001;...
                                                          3579     21631     36000     40857     52908     55270     74539    136286    150613    168000    188000    192137    216000    256450    278911,...
                                                        319800    321780    336000    372000    380000]';
info.R1151E.artfctdef.visual.R1151E_FR1_1_15Mar16_1654 = [1783543   2878942   2893342;...
                                                          1783805   2879084   2893577]';
info.R1151E.artfctdef.visual.R1151E_FR1_2_16Mar16_2141 = [377751   384001   594283     1492001;...     
                                                          381005   388000   595200     2543508]'; % 2nd half is messed up
info.R1166D.artfctdef.visual.R1166D_FR1_0_20Apr16_1846 = [    1    20126   116001   312549   331565   467468   502831   518791   564001   594416   607460   610940   620763   667472   716896   780001   844001   889682  1063750  1160396  1171113  1176364  1244001  1276920  1284001  1309521  1313480  1327028  1331903  1340001  1345642  1394097  1412485  1548755  1557190  1629343  1653271  1657472  1771424  1844001  1878706  1964670  1972001  2156001  2190601  2194674  2294484  2310537  2452001  2482863  2672303  2710706  2812815  2917073  2921706  3056001  3095730  3125726  3132001  3136001  3235383  3242371  3254593;...
                                                          16000    22759   117070   314570   333925   472881   512000   529731   589619   601473   607703   614788   628000   670445   723042   793002   846026   892000  1066413  1162788  1175638  1177518  1246163  1278509  1285481  1311094  1315082  1328603  1337671  1342614  1372000  1398266  1431412  1553066  1558509  1630868  1656000  1662913  1776413  1848000  1882933  1968000  1979046  2158155  2191932  2200000  2298526  2312000  2452679  2506622  2676000  2711529  2814499  2920000  2924000  3060000  3096141  3128000  3135050  3137574  3240000  3245651  3257143]';
info.R1166D.artfctdef.visual.R1166D_FR1_1_22Apr16_2012 = [21267    30359    56658    74464   194932   288477   511190   553037   740469   883686   993501  1093327  1194230  1220533  1224481  1300719  1304711  1496501  1506254  1654117  2126581  2162097  2268372  2364001  2465186  2469944  2485839  2576372  2692864  2793779;...
                                                          21586    31066    58659    80000   196000   289917   512558   557090   742316   884482   996445  1104000  1200000  1223537  1227042  1301893  1306018  1497312  1508603  1656312  2128000  2167598  2269590  2370806  2467296  2472974  2487356  2579138  2694066  2800986]'; 
info.R1166D.artfctdef.visual.R1166D_FR1_2_23Apr16_1248 = [   1    87621   498799   568001   592178   632001   732001   852061   874807   908834   910113  1040215  1112577  1233634  1244001  1253839  1291440  1304791  1430307  1603436  1616001  1672001  1704001  1716001  1760001  1926629  2033214  2052001  2076602  2116001  2280501  2405521  2532964  2664565  2679162  2689335  2697190  2784001  2889351  2894424  2948116  2949472  2952001  2956001  2990065  3062617  3124001  3155387  3170621;...
                                                          4000    92000   500413   572000   596000   637369   736000   854997   884000   908890   912857  1043231  1114651  1237514  1248000  1254110  1293739  1307130  1446368  1611435  1618759  1676000  1712000  1719392  1796000  1927912  2040000  2068000  2108000  2125393  2284000  2408000  2552000  2668000  2682465  2689671  2697453  2785437  2890987  2900000  2948660  2949711  2953953  2960000  2993058  3067267  3143243  3164000  3172125]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
                   
   
         


   
   
%%%% OLD 
%%% Sessions with multiple datasets in them
% R1124J_1_FR1_1_12Jan16_2013
% R1124J_1_FR1_1_12Jan16_2045
% 
% R1172E_FR1_0_09May16_1646
% R1172E_FR1_0_09May16_1754
% R1172E_FR1_0_09May16_1812
% R1172E_FR1_0_09May16_1827
% 
% R1164E_FR1_0_24Apr16_1643
% R1164E_FR1_0_24Apr16_1645
% 
% R1156D_FR1_3_21May16_1625 *
% R1156D_FR1_3_21May16_1633
% 
% R1150J_FR1_1_11Mar16_1606
% R1150J_FR1_1_11Mar16_1646
% 
% R1091N_PAL2_0a_11Jun15_1459
% R1091N_PAL2_0b_11Jun15_1612
% R1091N_PAL2_0c_11Jun15_1643


   

  














 



