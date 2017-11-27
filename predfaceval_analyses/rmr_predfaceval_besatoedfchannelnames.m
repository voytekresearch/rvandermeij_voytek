function label = rmr_predfaceval_besatoedfchannelnames(label)


%
% Add 'POL 's
for ichan = 1:numel(label)
  label{ichan} = ['POL ' label{ichan}];
end

%
% Add '-Ref's
for ichan = 1:numel(label)
  label{ichan} = [label{ichan} '-Ref'];
end

% Reverse renaming the annotions channel
ind = find(strcmp(label,'EDFann'));
if ~isempty(ind)
  label{ind} = 'EDF Annotations';
end

























