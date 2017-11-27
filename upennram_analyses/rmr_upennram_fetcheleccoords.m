function coords = rmr_upennram_fetcheleccoords(currsubj,label)

% fetch info
info = rmr_upennram_info;

% number of locmonts check
nlocmont = numel(info.sessions.(currsubj));
if nlocmont~=1
  error('more than one locmont?')
end
% localization/montage are currently always the first (zero-indexed)
loc  = '0';
mont = '0';


% fetch contacts from json
contactsjson = [info.datapath 'session_data/experiment_data/protocols/r1/subjects/' currsubj '/localizations/' loc '/montages/' mont '/neuroradiology/current_processed/contacts.json'];
if exist(contactsjson,'file')
  contacts = loadjson(contactsjson);
  contacts = contacts.(currsubj).contacts;
else
  info.montages.(currsubj){iuniloc} = NaN;
  warning(['contacts.json not found for ' currsubj ': ' contactsjson])
end

% extract MNI coords from contacts based on label
nchan  = numel(label);
coords = NaN(nchan,3);
for ichan = 1:nchan
  currcont = contacts.(label{ichan}).atlases.mni;
  coords(ichan,:) = [currcont.x currcont.y currcont.z];
  % safety check
  if ~isequal(contacts.(label{ichan}).code,label{ichan})
    error('seriously? goddamnit')
  end
end


