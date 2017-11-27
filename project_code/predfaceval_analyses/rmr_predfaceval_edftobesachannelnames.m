function label = rmr_predfacedval_edftobesachannelnames(label)


% 
% Remove all the 'POL 's
for ichan = 1:numel(label)
  if numel(label{ichan})>4 && strcmp(label{ichan}(1:4),'POL ')
    label{ichan} = label{ichan}(5:end);
  end
end

%
% Remove all the '-Ref's

for ichan = 1:numel(label)
  if numel(label{ichan})>4 && strcmp(label{ichan}(end-3:end),'-Ref')
    label{ichan} = label{ichan}(1:end-4);
  end
end

% shorthen the annotions channel
ind = find(strcmp(label,'EDF Annotations'));
if ~isempty(ind)
  label{ind} = 'EDFann';
end





























