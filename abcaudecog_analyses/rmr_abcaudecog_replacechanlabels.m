function data = rmr_abcaudecog_replacechanlabels(cfg,data)

%
%
% Replace channel labels by grid/strip labels
%
% cfg.subject = string, subject name
%        data = fieldtrip-style raw data structure
%
%
if isfield(data,'dimord') || ~isfield(data,'trial')
    error('only raw data structure allowed')
end


% fetch labels
label = data.label;
nchan = numel(data.label);

% replace label
unchanged = 0;
switch cfg.subject
    case {''}
        error('labels not defined yet')
        
    case 'GP15'
        % grid 1:64 = gdat 1:64
        % base LLATGD1:64
        for ichan = 1:nchan
            currlab = label{ichan}(5:end);
            switch str2num(currlab)
                case num2cell(1:64)
                    currlab = ['LLATGD' num2str(str2num(currlab),'%02d')]; % use '%02d' as option to add leading zero
                    label{ichan} = currlab;
                otherwise
                    unchanged = unchanged + 1;
            end
        end
        
    case 'GP22'
        % this ones relatively easy
        % grid 1:64 = gdat 17:32, 1:16, 33:64
        % base RLATGD1:64
        for ichan = 1:nchan
            currlab = label{ichan}(5:end);
            switch str2num(currlab)
                case num2cell(17:32)
                    currlab = ['RLATGD' num2str(str2num(currlab)-16,'%02d')]; % use '%02d' as option to add leading zero
                    label{ichan} = currlab;
                case num2cell(1:16)
                    currlab = ['RLATGD' num2str(str2num(currlab)+16)];
                    label{ichan} = currlab;
                case num2cell(33:64)
                    currlab = ['RLATGD' currlab];
                    label{ichan} = currlab;
                otherwise
                    unchanged = unchanged + 1;
            end
        end
        
    case 'GP28'
        % this one's easy
        % grid 1:64 = gdat 1:64
        % base RLATGD1:64
        for ichan = 1:nchan
            currlab = label{ichan}(5:end);
            if (str2num(currlab)<=64) && (str2num(currlab)>=1)
                currlab = ['RLATGD' num2str(str2num(currlab),'%02d')]; % use '%02d' as option to add leading zero
                label{ichan} = currlab;
            else
                unchanged = unchanged + 1;
            end
        end
        
    case 'GP35'
        % gdat     grid
        % 1-16     49-64
        % 17-48    17-48
        % 49-64     1-16
        % base LLATGD1:64
        for ichan = 1:nchan
            currlab = label{ichan}(5:end);
            switch str2num(currlab)
                case num2cell(1:16)
                    currlab = ['LLATGD' num2str(str2num(currlab)+48,'%02d')]; % use '%02d' as option to add leading zero
                    label{ichan} = currlab;
                case num2cell(17:48)
                    currlab = ['LLATGD' num2str(str2num(currlab),'%02d')];
                    label{ichan} = currlab;
                case num2cell(49:64)
                    currlab = ['LLATGD' num2str(str2num(currlab)-48,'%02d')];
                    label{ichan} = currlab;
                otherwise
                    unchanged = unchanged + 1;
            end
        end
        
    case 'ST1'
        % 1:64  = 1:64 RLATGD
        % 65:68 = 1:4 RIHST
        % 69:72 = 1:4 LIHST
        for ichan = 1:nchan
            currlab = label{ichan}(5:end);
            switch str2num(currlab)
                case num2cell(1:64)
                    currlab = ['RLATGD' num2str(str2num(currlab),'%02d')]; % use '%02d' as option to add leading zero
                    label{ichan} = currlab;
                case num2cell(65:68)
                    currlab = ['RIHST' num2str(str2num(currlab)-64,'%02d')];
                    label{ichan} = currlab;
                case num2cell(69:72)
                    currlab = ['LIHST' num2str(str2num(currlab)-68,'%02d')];
                    label{ichan} = currlab;
                otherwise
                    unchanged = unchanged + 1;
            end
        end
        
    case 'ST6'
        %  1-10 : LINAGD - 1-10 - LGA 1-10
        % 10-20 : RINAGD - 1-10 - RGA right 11-20
        % 21-30 : RINBGD - 1-10 - RGP right 1-10
        % 31-40 : LINBGD - 1-10 - LGP  11-20
        % 41-44 : RINTST - 1-4 - RIHS right 1-4
        % 45-48 : LINTST - 1-4 - LIHS 5-8
        % 49-56 : LANTST - 1-8 - LSA Anterior 1-8
        % 57-64 : LPOSST - 1-8 - LSP Posterior 1-8
        % 65-76 : LMEDST - 1-12 - LSM 1?12
        for ichan = 1:nchan
            currlab = label{ichan}(5:end);
            switch str2num(currlab)
                case num2cell(1:10)
                    currlab = ['LINAGD' num2str(str2num(currlab),'%02d')]; % use '%02d' as option to add leading zero
                    label{ichan} = currlab;
                case num2cell(11:20)
                    currlab = ['RINAGD' num2str(str2num(currlab)-10,'%02d')];
                    label{ichan} = currlab;
                case num2cell(21:30)
                    currlab = ['RINBGD' num2str(str2num(currlab)-20,'%02d')];
                    label{ichan} = currlab;
                case num2cell(31:40)
                    currlab = ['LINBGD' num2str(str2num(currlab)-30,'%02d')];
                    label{ichan} = currlab;
                case num2cell(41:44)
                    currlab = ['RINTST' num2str(str2num(currlab)-40,'%02d')];
                    label{ichan} = currlab;
                case num2cell(45:48)
                    currlab = ['LINTST' num2str(str2num(currlab)-44,'%02d')];
                    label{ichan} = currlab;
                case num2cell(49:56)
                    currlab = ['LANTST' num2str(str2num(currlab)-48,'%02d')];
                    label{ichan} = currlab;
                case num2cell(57:64)
                    currlab = ['LPOSST' num2str(str2num(currlab)-56,'%02d')];
                    label{ichan} = currlab;
                case num2cell(65:76)
                    currlab = ['LMEDST' num2str(str2num(currlab)-64,'%02d')];
                    label{ichan} = currlab;
                otherwise
                    unchanged = unchanged + 1;
            end
        end
        
        
    case 'ST8'
        %  1-16 = LANTGD 1-16  (left anterior temporal grid)
        % 17-32 = LPOSGD 17-32 (left posterior  temporal grid)
        % 33-42 = LAMYDP 1-10  (left amygdala depth)
        % 43-52 = LHIPDP 1-10  (left hippocampus depth)
        % 53-60 = LANTST 1-8   (left anterior strip)
        for ichan = 1:nchan
            currlab = label{ichan}(5:end);
            switch str2num(currlab)
                case num2cell(1:16)
                    currlab = ['LANTGD' num2str(str2num(currlab),'%02d')]; % use '%02d' as option to add leading zero
                    label{ichan} = currlab;
                case num2cell(17:32)
                    currlab = ['LPOSGD' num2str(str2num(currlab)-16,'%02d')];
                    label{ichan} = currlab;
                case num2cell(33:42)
                    currlab = ['LAMYDP' num2str(str2num(currlab)-32,'%02d')];
                    label{ichan} = currlab;
                case num2cell(43:52)
                    currlab = ['LHIPDP' num2str(str2num(currlab)-42,'%02d')];
                    label{ichan} = currlab;
                case num2cell(53:60)
                    currlab = ['LANTST' num2str(str2num(currlab)-52,'%02d')];
                    label{ichan} = currlab;
                otherwise
                    unchanged = unchanged + 1;
            end
        end
        
        
        
    case {'JH2','JH13'}
        % also kind of easy, channel names are in original filenames
        info = rmr_abcaudecog_info;
        chanfiles  = dir([info.datapath cfg.subject '/data/' info.session.(cfg.subject){1} '/' 'gdat*']);
        chanfiles  = {chanfiles.name};
        chanfiles  = chanfiles(:);
        nchanfiles = numel(chanfiles);
        % extract channel number and turn chanfiles into newlabels
        channum   = NaN(size(chanfiles));
        newlabel  = cell(size(chanfiles));
        for ichan = 1:nchanfiles
            currchan = chanfiles{ichan};
            prepos = strfind(currchan,'_');
            undpos = strfind(currchan,'_');
            prepos  = undpos(1);
            postpos = undpos(2);
            channum(ichan) = str2num(currchan((prepos+1):(postpos-1)));
            currnewlab = currchan(postpos+1:end-4);
            % in case of (S)EEG prefix, cut it off
            whisppos   = strfind(currnewlab,' ');
            if ~isempty(whisppos)
                currnewlab = currnewlab(whisppos+1:end);
            end
            % add zero to channel numbers below 10;
            if isempty(str2num(currnewlab(end-1:end)))
                currnewlab = [currnewlab(1:end-1) num2str(str2num(currnewlab(end)),'%02d')];
            end
            newlabel{ichan} = currnewlab;
        end
        % rename channels to be more informative
        for ichan = 1:nchanfiles
            currlab = newlabel{ichan};
            currprefix = currlab(1:end-2); % Only works when all labels have two number format (enforced above)
            switch currprefix
                % JH2 and JH13 mixed together
                case 'LH'
                    newlabel{ichan} = ['LLATGD' currlab(end-1:end)];
                case 'AD'
                    newlabel{ichan} = ['LANTDP' currlab(end-1:end)];
                case 'HD'
                    newlabel{ichan} = ['LHIPDP' currlab(end-1:end)];
                case {'LF','LFG'}
                    newlabel{ichan} = ['LFROGD' currlab(end-1:end)];
                case 'LMT'
                    newlabel{ichan} = ['LMEDST' currlab(end-1:end)];
                case 'LAT'
                    newlabel{ichan} = ['LANTST' currlab(end-1:end)];
                case 'LPT'
                    newlabel{ichan} = ['LPOSST' currlab(end-1:end)];
                case 'LOF'
                    newlabel{ichan} = ['LFROST' currlab(end-1:end)];
                case 'HDA'
                    newlabel{ichan} = ['LHIADP' currlab(end-1:end)];
                case 'HDB'
                    newlabel{ichan} = ['LHIBDP' currlab(end-1:end)];
                case 'IDA'
                    newlabel{ichan} = ['LINADP' currlab(end-1:end)];
                case 'IDB'
                    newlabel{ichan} = ['LINBDP' currlab(end-1:end)];
                case 'EKG'
                    % do nothing
                otherwise
                    error('woopsie')
            end
        end
        % replace old labels by new labels
        for ichan = 1:nchan
            currlab = label{ichan};
            % extract channel number
            currlabnum = str2num(currlab(5:end));
            % get index into newlabel
            [dum ind] = intersect(channum,currlabnum);
            if ~isempty(ind)
                label{ichan} = newlabel{ind};
            else
                unchanged = unchanged + 1;
            end
        end
        
    otherwise
        error('woopsie')
end


% resort labels and trials and put labels back in output
[label sortind] = sort(label);
data.label = label;
for itrial = 1:numel(data.trial)
    data.trial{itrial} = data.trial{itrial}(sortind,:);
end

% give warning when not all channels labels were changed
if unchanged>0
    warning([num2str(unchanged) ' channels did not have a mapping'])
end

















