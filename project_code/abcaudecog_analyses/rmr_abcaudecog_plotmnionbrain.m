function rmr_abcaudecog_plotmnionbrain

% get subject details
info = rmr_abcaudecog_info;

% discard subjects whose electrodes haven't been localized yet
[dum remind] = intersect(info.subj,{'GP14','GP15','ST8'});
info.subj(remind) = [];

% set n's, loop, plot locations
nsubj = numel(info.subj);
for isubj = 1:nsubj
    % set 
    currsubj = info.subj{isubj};
    
    % fetch locations
    mnicoords = rmr_abcaudecog_getmniloc(currsubj);
    
    % plot on mni brain
    figure('numbertitle','off','name',currsubj,'position',[ 1          49        1920         944])
    cfg = [];
    cfg.template     = 'surface_white_both';
    cfg.patchreduce  = 0.05;
    cfg.coordinates  = mnicoords;
    cfg.markersize   = 7;
    cfg.renderbrain  = 'no';
    roe_brainplot_chancircle(cfg);
end
