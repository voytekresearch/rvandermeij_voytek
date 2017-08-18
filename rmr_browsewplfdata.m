function rmr_browsewplfdata(cfg,wplfdata)

% RMR_BROWSEWPLFDATA allows for spatially browsing the results of Phase-Amplitude Coupling (PAC) analyses.
% Input data should be in the format of the output of rmr_phaseamplitudecoupling. In order to browse PAC,
% a FieldTrip-style electrode layout should also be provided. If not, an 'ordered' will be automatically
% generated, but this is not advised.
%
% When browsing, the electrode that is the seed is higlighted by a white rectangle. From a dropdown menu
% on the right it can be selected whether the seed is an amplitude- or phase-providing electrode. Wether
% PAC magnitude or PAC preferred phase of coupling is shown can be selected from the second dropdown box.
% To highlight the same frequencies in every electrode, you can use the highlight toggle. When clicked,
% draw a box on a electrode, and the same box will be visible on the others. This is especially usefully when
% visualizing PAC phase, as the phase for frequency-pairs outside the box will invisible.
% When viewing PAC magnitude the maximum of the color scale can be set manually or automatically. The minimum
% is always 0.
% In the bottom right a grand average is shown, from the currently selected seed only. This is mainly intended
% to visualize which frequency-pairs are which (the 'highlight' box will be visible here too).
%
% More information on the FieldTrip layout-structure can be found here: www.fieldtriptoolbox.org/tutorial/layout
% This also describes how to interactively create one using an image of the electrode grid.
% The electrode layout can also be generated semi-automatically using anatomical information using ft_prepare_layout.m.
% For this, a (cfg.)viewpoint for the brain should be selected, and only those electrodes that are visibile from this viewpoint.
% If anatomy is included using cfg.headshape and cfg.mri a brain outline will be included in the layout.
% Note, appending layouts of different electrode sets is current not yet supported by ft_appendlayout.
%
%
% Use as
%   rmr_browsewplfdata(cfg, wplfdata)
%
% The input data should be organised in a structure as obtained from RMR_PHASEAMPLITUDECOUPLING.
% The configuration structure contains options as specified below.
%
%        Input options:
%    cfg.layout      = a FieldTrip-style layout structure (if empty, an 'ordered' layout will be generated)
%    cfg.plotlabel   = 'yes' or 'no', plot the electrode labels on top of every PAC image.
%    cfg.plotoutline = 'yes' or 'no', whether or not to plot the outline as given in the layout structure
%    cfg.windowname  = string, title of figure window
%
%
% This function depends on the FieldTrip toolbox, which can be found at www.fieldtriptoolbox.org.
%
% Copyright (C) 2016-present, Roemer van der Meij
%

% Set defaults
cfg.layout       = ft_getopt(cfg, 'layout',       []);
cfg.windowname   = ft_getopt(cfg, 'windowname',   []);
cfg.plotlabel    = ft_getopt(cfg, 'plotlabel',    false);
cfg.plotoutline  = ft_getopt(cfg, 'plotoutline',  true);


% check whether a trial dimension is present in the data
if strcmp(wplfdata.dimord,'rpt_ampchan_phaschan_ampfreq_phasfreq')
  wplfdata.wplfrpt = wplfdata.wplf;
  wplfdata.wplf = squeeze(mean(wplfdata.wplfrpt,1));
  warning('averaging over repetitions')
end

%%% backwards compatability
if isfield(wplfdata,'wphaselockspctrm')
  wplfdata.wplf  = wplfdata.wphaselockspctrm;
  wplfdata.label = wplfdata.label_old;
  wplfdata = rmfield(wplfdata,{'wphaselockspctrm','label_old'});
  if ~isfield(wplfdata,'ampfreq') && ~isfield(wplfdata,'phasfreq')
    ampfreq  = wplfdata.freq_old(:);
    phasfreq = wplfdata.freq_old(:);
    wplfdata = rmfield(wplfdata,{'freq_old'});
  end
end
if strcmp(wplfdata.dimord,'part_ampchan_phaschan_ampfreq_phasfreq')
  if any(size(wplfdata.wplf) == 1)
    tmpsize = size(wplfdata.wplf);
    tmpsize = tmpsize(2:end);
    tmp = squeeze(mean(wplfdata.wplf,1));
    wplfdata.wplf = reshape(tmp,tmpsize);
    warning('averaging over random paritions (note: if normalization is over trials, this not the same as calculating wPLFs over all trials)')
  else
    wplfdata.wplf = squeeze(mean(wplfdata.wplf,1));
    warning('averaging over random paritions (note: if normalization is over trials, this not the same as calculating wPLFs over all trials)')
  end
end
%%% backwards compatability


% Set dimension sizes and other variables
nchan      = size(wplfdata.wplf,1);
label      = wplfdata.label;
ampfreq    = wplfdata.ampfreq(:);
phasfreq   = wplfdata.phasfreq(:);
nampfreq   = numel(ampfreq);
nphasfreq  = numel(phasfreq);
wplf       = wplfdata.wplf;


% prepare ticks and labels for the grand average
onefreqflg = false;
if nampfreq ~= 1
  amptickind   = 1:floor(nampfreq/6):nampfreq;
  ampticklabel = num2str(round(ampfreq(amptickind)));
else
  ampticklabel = num2str(round(ampfreq));
  amptickind   = 1;
  onefreqflg = true;
end
if nphasfreq ~= 1
  phastickind   = 1:floor(nphasfreq/6):nphasfreq;
  phasticklabel = num2str(round(phasfreq(phastickind)));
else
  phasticklabel = num2str(round(phasfreq));
  phastickind   = 1;
  onefreqflg = true;
end

% handle layout
if ~isempty(cfg.layout)
  lay = cfg.layout;
else
  % generate ordered layout
  cfglay = [];
  cfglay.layout = 'ordered';
  lay = ft_prepare_layout(cfglay,wplfdata);
end


% convert and scale layout to a 'stretched' format using [left bottom width height]
poslay = zeros(size(lay.pos,1),4);
outline = lay.outline;
roffset = 0.2;
loffset = 0.05;
boffset = 0.05;
toffset = 0.05;
% first insert the unscaled x,y coordinates
poslay(:,1:2) = lay.pos;
% % first, decenter x and y coords
% poslay(:,1) = poslay(:,1) -(lay.width/2);
% poslay(:,2) = poslay(:,2) -(lay.height/2);
% shift x coordinates to represent 'left offset'
for ioutline = 1:numel(outline)
  outline{ioutline}(:,1) = outline{ioutline}(:,1) - min(poslay(:,1));
end
poslay(:,1) = poslay(:,1) - min(poslay(:,1));
% shift y coordinates to represent 'bottom offset'
for ioutline = 1:numel(outline)
  outline{ioutline}(:,2) = outline{ioutline}(:,2) - min(poslay(:,2));
end
poslay(:,2) = poslay(:,2) - min(poslay(:,2));
% add width and height
poslay(:,3) = lay.width;
poslay(:,4) = lay.height;
% make em square by 0.9 times the width
poslay(:,3) = lay.width .* .9;
poslay(:,4) = lay.width .* .9;
% get scaling factor and normalize
scalefac = max(max(poslay(:,1:2) + poslay(:,3:4)));
for ioutline = 1:numel(outline)
  outline{ioutline} = outline{ioutline} ./ scalefac;
end
poslay = poslay ./ scalefac;
% add left offset
for ioutline = 1:numel(outline)
  outline{ioutline} = outline{ioutline} + loffset;
end
poslay(:,1:2) = poslay(:,1:2) + loffset;
% scale it such that the right side doesn't come close to 'offset' amount
scalefac = max(poslay(:,1) + poslay(:,3) - (1-roffset));
if scalefac > 0
  for ioutline = 1:numel(outline)
    outline{ioutline} = outline{ioutline} .* (1-scalefac);
  end
  poslay = poslay .* (1-scalefac);
end
% shift all to centre
distfromtop = (1-toffset) - max(poslay(:,2));
distfrombot = abs(min(poslay(:,2)) - boffset);
for ioutline = 1:numel(outline)
  outline{ioutline}(:,2) = outline{ioutline}(:,2) + ((distfromtop + distfrombot)/2);
end
poslay(:,2) = poslay(:,2) + ((distfromtop + distfrombot)/2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figure and set up guidat
if isempty(cfg.windowname)
  h = figure('name','wplfdata','numbertitle','off','units','normalized','position',[0 0.0208 1 0.8933]); % from a full screen fig window
else
  h = figure('name',cfg.windowname,'numbertitle','off','units','normalized','position',[0 0.0208 1 0.8933]); % from a full screen fig window
end
guidat = [];
guidat.wplf          = wplf;
guidat.lay           = lay;
guidat.poslay        = poslay;
guidat.seed          = round(nchan/2);
guidat.seeddim       = 2; % seed is amplitude (1) or phase(2) providing
guidat.wplfinfo      = 1; % display magnitude (1) or phase (2)
guidat.label         = label;
guidat.ampfreq       = ampfreq;
guidat.phasfreq      = phasfreq;
guidat.nchan         = nchan;
guidat.nampfreq      = nampfreq;
guidat.nphasfreq     = nphasfreq;
guidat.ampticklabel  = ampticklabel;
guidat.amptickind    = amptickind;
guidat.phasticklabel = phasticklabel;
guidat.phastickind   = phastickind;
guidat.onefreqflg    = onefreqflg;
guidat.magscalemax   = 0; % will cause auto determination on first draw
guidat.magscaledit   = 'off';
guidat.plotlabel     = istrue(cfg.plotlabel);
guidat.highlight     = [];
guidat.firstdraw     = 1;

% add guidata to figure
setappdata(h, 'guidat', guidat);

% plot outlines
if istrue(cfg.plotoutline)
  axis manual
  for ioutline = 1:numel(outline)
    xline = outline{ioutline}(:,1);
    yline = outline{ioutline}(:,2);
    line(xline,yline,'linewidth',3)
  end
  axis off
end

%%% set up ui functions
% set channel selection by click
set(h, 'WindowButtonupFcn', @select_seedchan_cb);

% set seed dimension switch
uicontrol('tag', 'seedsel', 'parent', h, 'units', 'normalized', 'style', 'popupmenu', 'string','Amplitude-prov. electrode|Phase-prov. electrode','value',guidat.seeddim,'position', [0.83 0.78 0.10 0.1],'callback',@select_seeddim_cb);
axes('position',[0.812 0.84 0.08 0.1],'tag', 'seedseltext1')
text(0.25,0.5,'Seed is:','tag','seedseltext1')
axis off
axes('position',[0.812 0.795 0.08 0.1],'tag', 'seedseltext2','parent',h)
text(0.25,0.5,'- seed provides phase frequencies','tag','seedseltext2')
axis off

% set wplf info switch
uicontrol('tag', 'wplfinfosel', 'parent', h, 'units', 'normalized', 'style', 'popupmenu', 'string','wPLF magnitude|wPLF phase','value',guidat.wplfinfo,'position', [0.83 0.7 0.10 0.1],'callback',@wplfinfosel_cb);
axes('position',[0.812 0.76 0.08 0.1],'tag', 'wplfinfoseltext')
text(0.25,0.5,'Showing:','tag','wplfinfoseltext')
axis off

% set editable wPLF magnitude scaling
uicontrol('tag', 'magsel', 'parent', h, 'units', 'normalized', 'style', 'edit', 'string','<N/A>','enable',guidat.magscaledit,'position', [0.83 0.72 0.035 0.02],'callback',@magsel_cb);
axes('position',[0.812 0.7025 0.08 0.1],'tag', 'magseltext')
text(0.25,0.5,'Scaling (magnitude) max:','tag','magseltext')
axis off
uicontrol('tag', 'magselauto', 'parent', h, 'units', 'normalized', 'style', 'checkbox', 'string','auto','value',1,'position', [0.83 0.70 0.035 0.02],'backgroundcolor',[.8 .8 .8],'callback',@magselauto_cb);

% set highlight box selection
uicontrol('tag', 'highlighttoggle', 'parent', h, 'units', 'normalized', 'style', 'togglebutton', 'string','highlight','value',0,'position',[0.83 0.65 0.06 0.025],'callback',@highlight_toggle_cb)
%uicontrol('tag', 'highlighttoggleauto', 'parent', h, 'units', 'normalized', 'style', 'checkbox', 'string','scale based on highlight','value',0,'position', [0.83 0.63 0.1 0.02],'backgroundcolor',[.8 .8 .8],'callback',);

% draw fig
redraw_cb(h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function highlight_toggle_cb(htoggle,event)
% get parent
parenth = get(htoggle,'parent');
% fetch current state and activate selection if needed
togglestate = get(htoggle,'value');
if togglestate==0
  set(parenth, 'WindowButtonupFcn', @select_seedchan_cb);
  set(parenth, 'WindowButtonDownFcn', []);
  set(parenth, 'WindowButtonMotionFcn',[]);
  set(htoggle,'backgroundColor',[0.8 0.8 0.8])
  set(parenth,'pointer','arrow')
  uiresume(parenth)
  % retrieve appdata
  guidat = getappdata(parenth, 'guidat');
  guidat.highlight = [];
  % add it to figure
  setappdata(parenth, 'guidat', guidat);
elseif togglestate==1
  set(parenth, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', true, 'clear', false, 'callback', {@highlight_rangesel_cb, parenth}, 'event', 'WindowButtonDownFcn'});
  set(parenth, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', true, 'clear', false, 'callback', {@highlight_rangesel_cb, parenth}, 'event', 'WindowButtonUpFcn'});
  set(parenth, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', true, 'clear', false, 'callback', {@highlight_rangesel_cb, parenth}, 'event', 'WindowButtonMotionFcn'});
  set(htoggle,'backgroundColor','g')
end

% redraw fig
redraw_cb(parenth)
uiresume(parenth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function highlight_rangesel_cb(h,range,varargin)
% retrieve appdata
guidat = getappdata(h, 'guidat');
% round axes coordinates and apply limits
range = round(range);
if range(2)>guidat.nphasfreq, range(2) = range(2)-1; end
if range(4)>guidat.nampfreq,  range(4) = range(4)-1; end
if range(1)<1,                range(1) = 1;          end
if range(3)<1,                range(3) = 1;          end
% add to guidata and reset options to not-select
guidat.highlight = range;
setappdata(h, 'guidat', guidat);
set(h, 'WindowButtonupFcn',    @select_seedchan_cb);
set(h, 'WindowButtonDownFcn',  []);
set(h, 'WindowButtonMotionFcn',[]);
set(h,'pointer','arrow')

% redraw fig
redraw_cb(h)
uiresume(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function magsel_cb(hedit,event)
% get parent
parenth = get(hedit,'parent');
% get current number
magscalemax = str2double(get(hedit,'string'));
% retrieve appdata
guidat = getappdata(parenth, 'guidat');
% put in scaling
guidat.magscalemax = magscalemax;
% add it to figure
setappdata(parenth, 'guidat', guidat);

% redraw fig
redraw_cb(parenth)
uiresume(parenth)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function magselauto_cb(hbox,event)
% get parent
parenth = get(hbox,'parent');
% get current option
val = get(hbox,'Value');
% retrieve appdata
guidat = getappdata(parenth, 'guidat');
% set editbox to enable if val = 0
magselh = findobj(get(parenth,'children'),'tag', 'magsel');
if val == 0
  set(magselh,'enable','on')
  guidat.magscaledit = 'on';
elseif val == 1
  set(magselh,'enable','off')
  guidat.magscaledit = 'off';
end
% add it to figure
setappdata(parenth, 'guidat', guidat);

% resume
uiresume(parenth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_seedchan_cb(h,event)

% get 'tag' from current gca
tag = get(gca,'tag');
% if tag of clicked element is longer than 13 elements, it's a channel!
if numel(tag) > 13
  % retrieve appdata
  guidat = getappdata(h, 'guidat');
  % put in seed
  guidat.seed = str2double(tag(14:end));
  % add it to figure
  setappdata(h, 'guidat', guidat);
  
  % redraw fig
  redraw_cb(h)
end
uiresume(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_seeddim_cb(hmenu,event)
% get parent and guidat
parenth = get(hmenu,'parent');
guidat = getappdata(parenth, 'guidat');

% apply selected option
val = get(hmenu,'Value');
guidat.seeddim = val;
% change accompanying text
delete(findobj(parenth,'tag', 'seedseltext2'));
if val == 1
  axes('position',[0.812 0.795 0.08 0.1],'tag', 'seedseltext2','parent',parenth)
  text(0.25,0.5,'- seed provides amplitude frequencies','tag','seedseltext2')
  axis off
elseif val == 2
  axes('position',[0.812 0.795 0.08 0.1],'tag', 'seedseltext2','parent',parenth)
  text(0.25,0.5,'- seed provides phase frequencies','tag','seedseltext2')
  axis off
end
% add it to figure
setappdata(parenth, 'guidat', guidat);

% draw fig
redraw_cb(parenth)
uiresume(parenth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wplfinfosel_cb(hmenu,event)
parenth = get(hmenu,'parent');
guidat  = getappdata(parenth, 'guidat');

% apply selected option
val = get(hmenu,'Value');
guidat.wplfinfo = val;
% turn scaling sel off if phase is being plotted
magselh     = findobj(get(parenth,'children'),'tag', 'magsel');
magselautoh = findobj(get(parenth,'children'),'tag', 'magselauto');
if val == 2
  % save current states
  guidata.magscalauto = get(magselautoh,'enable');
  set(magselh,'enable','off','string','<N/A>')
  set(magselautoh,'enable','off')
elseif val == 1
  set(magselh,'enable',guidat.magscaledit,'string',num2str(guidat.magscalemax))
  set(magselautoh,'enable','on')
end
% add it to figure
setappdata(parenth, 'guidat', guidat);

% draw fig
redraw_cb(parenth)
uiresume(parenth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUBFUNCTION to draw figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXME: modularize and remove duplication
function redraw_cb(h, eventdata)
guidat = getappdata(h, 'guidat');

% remove old drawn elements that can be refreshed by selections
delete(findobj(h,'tag', 'drawelem'));

% drawn based on whether amp/phase freqs are singular or not
% amp and phase freqs both consists of more than one freq --> imagesc

%%%%%%%%%%
% select data
if ~guidat.onefreqflg %%% both amp/phase freqs consist of more than one freq
  switch guidat.seeddim
    case 1
      switch guidat.wplfinfo
        case 1
          plotdat  = abs(permute(guidat.wplf(guidat.seed,:,:,:),[2 3 4 1]));
        case 2
          plotdat  = angle(permute(guidat.wplf(guidat.seed,:,:,:),[2 3 4 1]));
      end
    case 2
      switch guidat.wplfinfo
        case 1
          plotdat  = abs(permute(guidat.wplf(:,guidat.seed,:,:),[1 3 4 2]));
        case 2
          plotdat  = angle(permute(guidat.wplf(:,guidat.seed,:,:),[1 3 4 2]));
      end
  end
  
else %%% one freq is singular
  
  % select all bispectra from perspective of seed
  switch guidat.seeddim
    case 1
      switch guidat.wplfinfo
        case 1
          plotdat  = abs(squeeze(guidat.wplf(guidat.seed,:,:,:,:)));
        case 2
          plotdat  = angle(squeeze(guidat.wplf(guidat.seed,:,:,:,:)));
      end
    case 2
      switch guidat.wplfinfo
        case 1
          plotdat  = abs(squeeze(guidat.wplf(:,guidat.seed,:,:,:)));
        case 2
          plotdat  = angle(squeeze(guidat.wplf(:,guidat.seed,:,:,:)));
      end
  end
end
%%%%%%%%%%


%%%%%%%%%%
% determine clim/ylim
if  ~guidat.onefreqflg  %%% both amp/phase freqs consist of more than one freq
  
  if guidat.wplfinfo == 1
    % determine if auto-scaling or not, and do it
    magselh = findobj(get(h,'children'),'tag', 'magsel');
    if strcmp(get(magselh,'enable'),'off')
      if all(all(all(plotdat>=0)))
        clim = [0 max(max(max(plotdat)))];
      else
        clim = [-max(max(max(abs(plotdat)))) max(max(max(abs(plotdat))))];
      end
      clim(1) = str2double(num2str(clim(1),'%.2f')); % round to two decimals
      clim(2) = str2double(num2str(clim(2),'%.2f')); % round to two decimals
      % update scaling
      set(magselh,'string',num2str(clim(2)))
      guidat.magscalemax = clim(2);
    else
      if all(all(all(plotdat>=0)))
        clim = [0 guidat.magscalemax];
      else
        clim = [-guidat.magscalemax guidat.magscalemax];
      end
    end
    colormap(parula)
  else
    clim = [-pi pi];
    colormap(hsv)
  end
  
else %%% one freq is singular
  
  if guidat.wplfinfo == 1
    % determine if auto-scaling or not, and do it
    magselh = findobj(get(h,'children'),'tag', 'magsel');
    if strcmp(get(magselh,'enable'),'off')
      ylim = [0 max(max(max(plotdat)))];
      ylim(2) = str2double(num2str(ylim(2),'%.2f')); % round to two decimals
      % update scaling
      set(magselh,'string',num2str(ylim(2)))
      guidat.magscalemax = ylim(2);
    else
      ylim = [0 guidat.magscalemax];
    end
  else
    ylim = [-pi pi];
  end
end
%%%%%%%%%%


%%%%%%%%%%
% plot each channel
for ichan = 1:guidat.nchan
  % find/set axes
  if guidat.firstdraw
    axes('position',guidat.poslay(ichan,:));
    set(gca,'tag',['drawelemichan' num2str(ichan)])
    hold on
  else
    currah = findobj(h,'tag',['drawelemichan' num2str(ichan)]);
    axes(currah);
  end
  
  % plot wplfs
  if  ~guidat.onefreqflg %%% both amp/phase freqs consist of more than one freq
    
    currplotdat = squeeze(plotdat(ichan,:,:));
    if ~isempty(guidat.highlight) && (guidat.wplfinfo == 2)
      % (think of it as the data having an exact range of min=clim(1) to max=(clim2), convert this range to 0-64) FIXME: This isn't indentical to transparancy, but it's close
      cdat = (currplotdat + -clim(1)) * (64 / (-clim(1) + clim(2)));
      range = guidat.highlight;
      satmask = zeros(size(currplotdat));
      satmask(range(3):range(4),range(1):range(2)) = 1;
      satmask = (satmask .* .9) + .1;
      % ind->rgb->hsv ||change saturation values||  hsv->rgb ->  plot
      rgbcdat = ind2rgb(uint8(floor(cdat)), colormap);
      hsvcdat = rgb2hsv(rgbcdat);
      hsvcdat(:,:,2) = hsvcdat(:,:,2) .* satmask;
      rgbcdatsat = hsv2rgb(hsvcdat);
      currplotdat = rgbcdatsat;
    end
    currih = imagesc(currplotdat,clim);
    set(currih,'tag','drawelem')
    axis tight
    axis xy
  
  else %%% one freq is singular
    
    % plot wplfs
    currph = plot(squeeze(plotdat(ichan,:,:)),'color',[.2 .6 1]);
    set(currph,'tag','drawelem')
    axis tight
    set(gca,'ylim',ylim)
    if guidat.nampfreq == 1
      view(0,90)
    elseif guidat.nphasfreq == 1
      view(90,-90)
    end
  end
  
  % remove ticks
  set(gca,'YTick', [], 'XTick', [])
  
  % plot label if requested
  if guidat.plotlabel && guidat.firstdraw
    title(guidat.lay.label{ichan})
  end
  
  % plot a box around a highlight if desured
  if ~isempty(guidat.highlight) && (guidat.wplfinfo == 1)
    range = guidat.highlight;
    xcoords = [range(1) range(2) range(2) range(1) range(1)];
    ycoords = [range(3) range(3) range(4) range(4) range(3)];
    line(xcoords, ycoords,'color',[.4 .4 .4], 'linewidth', 1,'tag','drawelem')
  end
end
%%%%%%%%%%


% plot a box around selected seed
seedpos = guidat.poslay(guidat.seed,:);
seedshift = 0.013;
seedpos = seedpos + [-seedshift -seedshift 2.*seedshift 2.*seedshift];
axes('position',seedpos)
axis off
line([0.1 .9 .9 .1 .1],[0.1 0.1 0.9 0.9 .1],'color',[.7 .7 .7],'linewidth',4,'tag','drawelem')
set(gca,'tag','drawelem')

%%%%%%%%%%
% place a grand average with axes ticks on the right
% plot based on layout
prop = guidat.poslay(1,4) ./ guidat.poslay(1,3);
currplotdat = angle(squeeze(mean(exp(1i*plotdat),1)));
axes('position',[0.83 0.3 0.15 0.15.*prop])
if ~guidat.onefreqflg %%% both amp/phase freqs consist of more than one freq  
  
  if ~isempty(guidat.highlight) && (guidat.wplfinfo == 2)
    % (think of it as the data having an exact range of min=clim(1) to max=(clim2), convert this range to 0-64) FIXME: This isn't indentical to transparancy, but it's close
    cdat = (currplotdat + -clim(1)) * (64 / (-clim(1) + clim(2)));
    range = guidat.highlight;
    satmask = zeros(size(currplotdat));
    satmask(range(3):range(4),range(1):range(2)) = 1;
    satmask = (satmask .* .9) + .1;
    % ind->rgb->hsv ||change saturation values||  hsv->rgb ->  plot
    rgbcdat = ind2rgb(uint8(floor(cdat)), colormap);
    hsvcdat = rgb2hsv(rgbcdat);
    hsvcdat(:,:,2) = hsvcdat(:,:,2) .* satmask;
    rgbcdatsat = hsv2rgb(hsvcdat);
    currplotdat = rgbcdatsat;
  end
  imagesc(currplotdat,clim)
  axis tight;
  axis xy
  set(gca,'YTick', guidat.amptickind ,'YTickLabel', guidat.ampticklabel,'XTick', guidat.phastickind,'XTickLabel', guidat.phasticklabel)
  xlabel('phase frequency (Hz)')
  ylabel('amplitude frequency (Hz)')
  set(gca,'tag','drawelem')
  caxis(clim);
  title('Grand average of selected seed')
  
  % place a colorbar
  axes('position',[0.83 0.1 0.15 0.15])
  caxis(clim);
  hc = colorbar('south','tag','drawelem');
  axis off
  if guidat.wplfinfo == 2
    set(hc,'ticks',[-pi 0 pi],'ticklabels',{'-pi','0','pi'})
    text(0.5,0.5,'wPLF phase','tag','drawelem','horizontalalignment','center')
  else
    text(0.5,0.5,'wPLF magnitude','tag','drawelem','horizontalalignment','center')
  end
  
else %%% one freq is singular
  
  plot(angle(squeeze(mean(exp(1i*plotdat),1))),'color',[.2 .6 1]);
  axis tight;
  set(gca,'ylim',ylim)
  if guidat.wplfinfo == 2
    set(gca,'ytick',[-pi 0 pi],'yticklabel',{'-pi','0','pi'})
    ylabel('wPLF phase');
  else
    ylabel('wPLF magnitude');
  end
  
  if guidat.nampfreq == 1
    set(gca,'XTick', guidat.phastickind,'XTickLabel', guidat.phasticklabel)
    xlabel('phase frequency (Hz)');
    view(0,90)
  elseif guidat.nphasfreq == 1
    set(gca,'XTick', guidat.amptickind,'XTickLabel', guidat.ampticklabel)
    xlabel('amplitude frequency (Hz)');
    view(90,-90)
  end
  set(gca,'tag','drawelem')
  title('Grand average selected seed')
    
  % mention scaling
  axes('position',[0.83 0.1 0.15 0.15])
  text(0.25,0.5,['Scaling of all wPLFs: ' num2str(ylim(1)) ' <-> ' num2str(round(ylim(2).*1000)./1000)],'tag','drawelem')
  axis off
end
%%%%%%%%%%


% plot a box around a highlight if present
if ~isempty(guidat.highlight) && (guidat.wplfinfo == 1)
  range = guidat.highlight;
  xcoords = [range(1) range(2) range(2) range(1) range(1)];
  ycoords = [range(3) range(3) range(4) range(4) range(3)];
  line(xcoords, ycoords,'color',[.5 .5 .5], 'linewidth', 1,'tag','drawelem')
end


% flg
guidat.firstdraw = 0;

% add it to figure
setappdata(h, 'guidat', guidat);

% resume
uiresume(h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








