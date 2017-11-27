function mnicoords = rmr_abcaudecog_getchaninfo(subj)

% 
% 
% Reads electrode info out of .mgrid files
% 
% 
% 

% read data details
info = rmr_abcaudecog_info;

% get location of mgrid files for each subject
fn = [info.datapath subj '/Recon/MNI_3D/'];
switch subj
    case {'GP14','GP15','ST8'}
        error('no localization yet')
    case 'GP22'
        fn = [fn 'trans_NLbet_Rigid_grid.mgrid'];
    case 'GP28'
        fn = [fn 'MNI_trans_grid.mgrid'];
    case 'GP35'
        fn = [fn 'MNI_trans_Pointreg_grid.mgrid'];
    case 'JH2'
        fn = [fn 'MNI_JH2_grid_MR.mgrid'];
    case 'JH13'
        fn = [fn 'JH13_trans_grids.mgrid'];
    case 'ST1'
        fn = [fn 'MNI_grid.mgrid'];
    case 'ST6'
        fn = [fn 'ST6_MNI.mgrid'];
    otherwise
        error('woopsie')
end


% open file and grab all text, put it in a cell-array, 1 cell per line of text
fileid = fopen(fn);
grab = textscan(fileid, '%s','delimiter','/n');
fclose(fileid);
grab = grab{1};
%grab = fscanf(fileid, '%10000c'); % this works, but it's cheating. %10000c should be a proper wildcard

% search grab for entries for position, fetch position, and parse
posind = find(strcmp('#Position',grab));
pos    = grab(posind+1);
for ipos = 1:numel(pos)
    pos{ipos} = str2num(pos{ipos});
end
mnicoords = cat(1,pos{:});

% perform sanity check by making sure we have coordinates for each electrode in the mgrid file
% fetch number of mentions of electrodes
elecind = find(strncmp('# Electrode',grab,11));
elecind = elecind(~strncmp('# Electrode Grid',grab(elecind),16));
if numel(elecind) ~= size(mnicoords,1)
    error(['number of electrodes for which mni coordinates were extracted is not equal to number of electrodes in mgrid file for subjec ' subj])
end

% apply transform to get locations actually in MNI space
switch subj
    %case {'GP14','GP15','ST8'}
    case {'ST6'}
        % cropped error
        mnicoords(:,1) = -(mnicoords(:,1)-77);
        mnicoords(:,2) = -(mnicoords(:,2)-73);
        mnicoords(:,3) = (mnicoords(:,3)-60);
    case {'GP22','GP28','GP35','JH2','JH13','ST1'}
        % standard error
        mnicoords(:,1) = -(mnicoords(:,1)-90);
        mnicoords(:,2) = -(mnicoords(:,2)-90);
        mnicoords(:,3) = (mnicoords(:,3)-72);
    otherwise
        error('woopsie')
end



















