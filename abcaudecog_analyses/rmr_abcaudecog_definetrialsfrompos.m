function trl = rmr_abcaudecog_definetrialsfrompos(cfg)

%
%
% Read event files (.pos) and
%
% Input:
%  cfg.subj    = subject name
%  cfg.isess   = session NUMBER
%  cfg.fsample = sampling rate, used to compute RT
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                                         %%%%%%%
%%%%%%  For event codes see event_codes.xls from Aurelie                       %%%%%%%
%%%%%%                                                                         %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Time between successive sounds is 400ms to 600ms.
% In the LR task (all monaural, both in high/low pitch, and L and R attented):
%    unattended standards
%    unattended deviants
%    attended standards
%    attended deviants (targets)
%    unattend binaurals
% In the binaural task
%    monaural standards
%    binaural standards
%
%
% This function finds all trials
%   1)  4th column of trl: condition (attL = 1, attR = 2, attB = 3)
%   2)  5th column of trl: standard/deviant (std = 1, dev = 2, 3rd cond = NaN)
%   3)  6th column of trl: target/no-target (1/0)
%   4)  7th column of trl: hit/miss/FA/CR (1/2/3/4)
%   5)  8th column of trl: RT
%   6)  9th column of trl: left/right/bin stim (left = 1, right = 2, binaural = 3)
%   7) 10th column of trl: low/high pitch (1/2)
%   8) 11th column of trl: session/block number
%   9) 12th column of trl: trial number in order of appearance
%  10) 13th column of trl: event code
%
%

if ~isfield(cfg,'subj') || ~isfield(cfg,'isess') || ~isfield(cfg,'fsample')
    error('improper input')
end


% get basics
info = rmr_abcaudecog_info;
currsess = lower(info.session.(cfg.subj){cfg.isess});
switch cfg.subj
    case {'JH2','JH13'}
        posfn = [info.datapath 'from_aurelie/pos_files/' lower(cfg.subj) '/'  lower(cfg.subj) '.asatt.pos'];
    case {'GP15'}
        undpos   = findstr(currsess,'_');
        currsess = currsess(undpos+1:end);
        posfn = [info.datapath 'from_aurelie/pos_files/' lower(cfg.subj) '/'  lower(cfg.subj) '.asatt.' currsess '.3052.pos'];
    case {'GP22','GP35'}
        undpos   = findstr(currsess,'_');
        currsess = currsess(undpos+1:end);
        posfn = [info.datapath 'from_aurelie/pos_files/' lower(cfg.subj) '/'  lower(cfg.subj) '.asatt.' currsess '.24414.cor.res3052.pos'];
    case {'GP28'}
        undpos   = findstr(currsess,'_');
        currsess = currsess(undpos+1:end);
        posfn = [info.datapath 'from_aurelie/pos_files/' lower(cfg.subj) '/'  lower(cfg.subj) '.asatt.' currsess '.24414.res3052.pos'];
    case {'ST1'}
        undpos   = findstr(currsess,'_');
        currsess = currsess(undpos+1:end);
        posfn = [info.datapath 'from_aurelie/pos_files/' lower(cfg.subj) '/'  lower(cfg.subj) '.asatt.' currsess '.3052.cor.pos'];
    case {'ST6'}
        undpos   = findstr(currsess,'_');
        currsess = ['b' currsess(undpos+1:end)];
        if strcmp(currsess,'b09')
            currsess = 'b9';
        end
        posfn = [info.datapath 'from_aurelie/pos_files/' lower(cfg.subj) '/'  lower(cfg.subj) '.asatt.' currsess '.24414.cor.res3052.pos'];
    case {'ST8'}
        zeropos  = findstr(currsess,'0');
        currsess = ['b' currsess(zeropos+1:end)];
        posfn = [info.datapath 'from_aurelie/pos_files/' lower(cfg.subj) '/'  lower(cfg.subj) '.asatt.' currsess '.24414.cor.res3052.pos'];
    otherwise
        error('woopsie!')
end


% open file and parse
fileid = fopen(posfn);
grab = textscan(fileid, '%f %f %f' ,'delimiter','/n');
fclose(fileid);
eventsample = grab{1};
eventcode  = grab{2};

% find trials and build trl
trl = [];
trlnum = 0;
multreptrl = 0;
allstimcodes = [111 113   112 114   221 223   222 224   211 213   212 214   121 123   122 124   11 13   12 14   21 23   22 24   15 25 135];
nevent = numel(eventcode);
for ievent = 1:nevent
    
    % select trials using all stimulus codes
    if ismember(eventcode(ievent),allstimcodes)
        
        % check whether at least two response were made, and discard if so
        if ievent <= nevent-2;
            if (eventcode(ievent+1) == 10) && (eventcode(ievent+2) == 10)
                % count
                multreptrl = multreptrl + 1;
                continue
            end
        end
        % discard trial if event occurred before it
        if ievent > 1
            if eventcode(ievent-1) == 10
                continue
            end
        end
        
        % check whether stimulus was a target, and wether a response was present afterwards
        istarget   = ismember(eventcode(ievent),[113 114 123 124 135]);
        if ievent < nevent
            hadresp = eventcode(ievent+1) == 10;
        else
            hadresp = 0;
        end
        
        % check whether stimulus was a standard or deviant
        if     ismember(eventcode(ievent),[111   112   221   222   211   212   121   122   11   12   21   22])
            standev = 1; % standard
        elseif ismember(eventcode(ievent),[113   114   223   224   213   214   123   124   13   14   23   24])
            standev = 2; % deviant
        elseif ismember(eventcode(ievent),[15 25 135])
            standev = NaN; % binaural sounds, can't be classified as standard or deviant
        else
            error('woopsie!')
        end
        
        
        % compute RT and set response status
        if hadresp
            RT =  round(((eventsample(ievent+1)-eventsample(ievent)) ./ cfg.fsample) * 1000); % in ms
        else
            RT = 0;
        end
        if istarget && hadresp
            respstat = 1; % HIT
        elseif istarget && ~hadresp
            respstat = 2; % MISS
        elseif ~istarget && hadresp
            respstat = 3; % FA
        else % ~istarget && ~hadresp
            respstat = 4; % CR
        end
        
        % set other trial classifiers
        % high low pitch
        if     ismember(eventcode(ievent),[111 113  221 223  211 213   121 123   11 13   21 23])
            lohipitch = 1; % low pitch
        elseif ismember(eventcode(ievent),[112 114   222 224   212 214   122 124   12 14   22 24])
            lohipitch = 2; % high pitch
        elseif ismember(eventcode(ievent),[15 25 135])
            lohipitch = NaN; % binaural sounds, can't be classified as low or high pitch
        else
            error('woopsie!')
        end
        
        % condition lrbin
        if     ismember(eventcode(ievent),[111 113   112 114   221 223   222 224   15])
            condlrbin = 1; % left attend condition
        elseif ismember(eventcode(ievent),[211 213   212 214   121 123   122 124   25])
            condlrbin = 2; % right attend condition
        elseif ismember(eventcode(ievent),[11 13   12 14   21 23   22 24  135])
            condlrbin = 3;  % binaural attend condition
        else
            error('woopsie!')
        end
        % stimulus lrbin
        if     ismember(eventcode(ievent),[111 113   112 114   211 213   212 214   11 13   12 14])
            stimlrbin = 1; % left stimulus
        elseif ismember(eventcode(ievent),[221 223   222 224   121 123   122 124   21 23   22 24])
            stimlrbin = 2; % right stimulus
        elseif ismember(eventcode(ievent),[15 25 135])
            stimlrbin = 3;  % binaural stimulus
        else
            error('woopsie!')
        end
        
        % find next stimulus event, and use it as the boundary of the current trial
        % if no next stimulus event, make the trial 600ms
        % if next event occurs after 600 ms, bound to 600ms (at least necessary for last/first stimulus of different blocks)
        foundnextstim = false;
        nextstimevent = ievent+1; % start from now
        while ~foundnextstim && nextstimevent <= nevent
            if ismember(eventcode(nextstimevent),allstimcodes)
                foundnextstim = true;
            else
                nextstimevent = nextstimevent + 1;
            end
        end
        % check for stimuli that occur to fast after the previous one (which happend in ST1, new stim after a few samples, and turned out GP35 had some bad ones too)
        % ISI should have lower bound of 300
        if foundnextstim
            if ((eventsample(nextstimevent)-eventsample(ievent)) ./ cfg.fsample) < 0.290 % use 290ms for some leeway, in case of timing inaccuracies
                continue
            end
        end
        % set samples
        begsample = eventsample(ievent);
        if foundnextstim
            if (eventsample(nextstimevent)-1) <= (begsample + floor(0.6 * cfg.fsample))
                endsample = eventsample(nextstimevent)-1;
            else
                endsample = begsample + floor(0.6 * cfg.fsample);
            end
        else
            endsample = begsample + floor(0.6 * cfg.fsample);
        end
        offset = 0;
        
        % build trl
        trlnum = trlnum + 1; % increment
        trl(end+1,:) = [begsample endsample offset condlrbin standev istarget respstat RT stimlrbin lohipitch cfg.isess trlnum eventcode(ievent)];
    end
end

% provide info
disp(['detected ' num2str(trlnum) ' trials, ' num2str(multreptrl) ' trials had more than one response and were removed'])













