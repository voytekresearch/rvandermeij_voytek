function rmr_upennram_getfourier





% fetch info
info = rmr_upennram_info;
info.subj = info.subjselmains;

%info.savepath = '/Users/roemer/Work/Data/tmpdata/upennram/';

for     isubj = 1  :numel(info.subj)
  % set
  currsubj   = info.subj{isubj};
  
  % for artifacts been detected
  if ~info.(currsubj).trlartfctflg
    disp(['artifacts not yet detected for ' currsubj])
    continue
  end
  
  
  % check whether already done
  fn = [];
  fn{1} = [info.savepath currsubj '_' 'fourier_3to40hz_3cyc'  '.mat'];
  fn{2} = [info.savepath currsubj '_' 'dataetc_3to40hz_3cyc'  '.mat'];
  if  all(cellfun(@exist,fn(:),repmat({'file'},[numel(fn) 1])))
    continue
  end
  
  % fetch and preprocess data
  cfg = [];
  cfg.useenc         = true;
  cfg.userec         = false;
  cfg.encduration    = 1.6;
  cfg.encprestim     = 0;
  cfg.encpoststim    = 0;
  cfg.demean      = 'yes';
  cfg.detrend     = 'yes';
  cfg.reref       = 'yes';
  cfg.bsfilter    = 'yes';
  cfg.lpfilter    = 'yes';
  data = rmr_upennram_fetchpreprocessdata(cfg,currsubj,info);
  
  
  
  % get fourier
  % get timeoi's for welch tapering
  freqoi    = 3:1:40;
  taplength = 3./freqoi;
  % get them
  [dum, maxind] = max(cellfun(@numel,data.time));
  timeoi = data.time{maxind};
  % get fourier
  fourier = rmr_SPACEinput_electrophysiology(data,freqoi,timeoi,taplength,'hanning',[]);
  
  % save fourier
  fourfn = fn{1};
  save(fourfn,'fourier','-v7.3')
  % save data and stuff
  data = rmfield(data,'trial');
  datafn = fn{2};
  save(datafn,'data','freqoi','timeoi')
end % exist






































