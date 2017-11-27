function roe_brainplot_chancircle(cfg)

% ROE_BRAINPLOT_CHANCIRCLE plots channels on a template brain
% by only requiring coordinates, and a viewpoint
% This function use the current figure window(!)
%
%
% Use as
%    roe_brainplot_chancircle(cfg)
%
%
% Necessary:
%    cfg.coordinates   = chanX3 matrix of X,Y,Z coordinates
%    cfg.viewpoint     = 1X2 vector, viewpoint of the brain
%                       viewpoint:   [90  0] - RSAG
%                                    [90  0] - RSAG_INT
%                                    [0   0] - Superior view
%                                    [-90 0] - LSAG
%                                    [-90 0] - LSAG_INT
%                                    [-180 -90] - Inferior view
%    cfg.displace      = 1X3 vector, X,Y,Z displacement of electrodes to place them on top of brain
%                       displacement: [ 100    0    0] - RSAG
%                                     [ 100    0    0] - RSAG_INT
%                                     [   0 -100    0] - Superior view
%                                     [-100    0    0] - LSAG
%                                     [-100    0    0] - LSAG_INT
%                                     [   0    0 -100] - Inferior view
%
%
%
% Optional:
%    cfg.alpha         = default is 1
%    cfg.braincolor    = color-vector
%    cfg.markersize    = default is 3 (can also be 1xnchan cell aray or 1xnchan vector containing sizes)
%    cfg.markercolor   = color-vector, color-string, or 'cparam' (can also be 1xnchan cell aray containing color vectors/strings)
%                        when using 'cparam', cfg.cparam also needs to be specified
%    cfg.markersymbol  = string containing marker (can also be 1xnchan cell aray containing markers)
%    cfg.cparam        = vector to use for coloring markers (overwrites cfg.markercolor)
%                        should have dimensions channels*1
%    cfg.clim          = to be used in conjuction with cfg.cparam, can be [cmin cmax], 'absmax' or 'minmax'
%                        (default is 'absmax')
%    cfg.colormap      = colormap to be used when plotting markers based on 'cparam'
%    cfg.colorbar      = 'yes' or 'no', displays colorbar
%    cfg.label         = 1xNchan cell-array or 'no'
%    cfg.labelsize     = fontsize
%    cfg.labelweigth   = font weight
%    cfg.labelcolor    = font color
%    cfg.patchreduce   = 'no', or a reduction-factor (e.g. 0.01)
%    cfg.template      = 'wholebrain.mat', 'leftbrain.mat', 'rightbrain.mat' (all TAL), 'surface_white_both' (MNI)
%

% Hidden
% cfg.axissqueeze = 'yes', or 'no'
% cfg.plotelec    = 'yes', or 'no'

% Set defaults
if ~isfield(cfg, 'alpha'),           cfg.alpha         = 1;                        end
if ~isfield(cfg, 'braincolor'),      cfg.braincolor    =  [.75 .75 .75];           end
if ~isfield(cfg, 'markersize'),      cfg.markersize    = 60;                       end
if ~isfield(cfg, 'markercolor'),     cfg.markercolor   = 'g';                      end
if ~isfield(cfg, 'markersymbol'),    cfg.markersymbol  = 'o';                      end
if ~isfield(cfg, 'clim'),            cfg.clim          = 'maxmin';                 end
if ~isfield(cfg, 'colormap'),        cfg.colormap      = 'jet';                    end
if ~isfield(cfg, 'colorbar'),        cfg.colorbar      = 'no';                     end
if ~isfield(cfg, 'label'),           cfg.label         = 'no';                     end
if ~isfield(cfg, 'labelsize'),       cfg.labelsize     = 8;                        end
if ~isfield(cfg, 'labelweigth'),     cfg.labelweigth   = 'normal';                 end
if ~isfield(cfg, 'labelcolor'),      cfg.labelcolor    =  [1 1 1];                 end
if ~isfield(cfg, 'patchreduce'),     cfg.patchreduce   = 'no';                     end
if ~isfield(cfg, 'viewpoint'),       cfg.viewpoint     = [90 0];                   end
if ~isfield(cfg, 'displace'),        cfg.displace      = [0 0 0];                  end
if ~isfield(cfg, 'coordinates'),       error('coordinates needed'),                end
if ~isfield(cfg, 'renderbrain'),     cfg.renderbrain   = 'yes';                    end
if ~isfield(cfg, 'template'),        cfg.template      = 'wholebrain.mat';         end
if ~isfield(cfg, 'axissqueeze'),     cfg.axissqueeze   = 'no';                     end
if ~isfield(cfg, 'plotelec'),        cfg.plotelec      = 'yes';                     end


% set nchan and coords
nchan = size(cfg.coordinates,1);
coords = cfg.coordinates;
coords = coords + repmat(cfg.displace,[size(coords,1) 1]);


% Set coloring of markers when using cparam
chancolor = [];
if isequal(cfg.markercolor,'cparam') && isfield(cfg,'cparam')
  % Get clim
  if strcmp(cfg.clim,'maxmin')
    cmin = min(cfg.cparam);
    cmax = max(cfg.cparam);
  elseif strcmp(cfg.clim,'absmax')
    cmin = -max(max(abs(cfg.cparam)));
    cmax = max(max(abs(cfg.cparam)));
  else
    cmin = cfg.clim(1);
    cmax = cfg.clim(2);
  end
  
  % Create colormat based on colormap (make it huge for indexing precision)
  colormat = eval([cfg.colormap '(1e5);']);
  
  % Transform cparam to go from 0 1 and multiply with 10e6 to act as index into colormat
  cparam = round(((cfg.cparam + -cmin) .* (1 / (-cmin + cmax))) .* 10e4);
  cparam(cparam==0) = 1; % if a color index is zero, zet it to the first color-value (which has index 1)
  cparam(cparam<0) = 1;
  cparam(cparam>1e5) = 1e5;
  colormat(1,:) = [.5 .5 .5];
  colormat(1e5,:) = [.5 .5 .5];
  for ichan = 1:nchan
    chancolor{ichan} = colormat(cparam(ichan),:);
  end
  cfg.markercolor = chancolor;
end

% Transform markercolor to cell-array when not using cparam
if ~iscell(cfg.markercolor) && (size(cfg.markercolor,1)==1) && ~isfield(cfg,'cparam')
  tempcolor = cfg.markercolor;
  cfg.markercolor = [];
  for ichan = 1:nchan
    cfg.markercolor{ichan} = tempcolor;
  end
elseif ~(size(cfg.markercolor,1)==1) && ~(size(cfg.markercolor,1)==nchan)
  error('cfg.markercolor should either have length = 1 or length = nchan')
end

% Transform markersymbol to cell-array
if ~iscell(cfg.markersymbol)
  tempsymbol = cfg.markersymbol;
  cfg.markersymbol = [];
  for ichan = 1:nchan
    cfg.markersymbol{ichan} = tempsymbol;
  end
elseif ~(length(cfg.markersymbol)==1) && ~(length(cfg.markersymbol)==nchan)
  error('cfg.markersymbol should either have length = 1 or length = nchan')
end

% Transform markersize to cell-array
if ~iscell(cfg.markersize)
  tempsize = cfg.markersize;
  cfg.markersize = [];
  for ichan = 1:nchan
    cfg.markersize{ichan} = tempsize;
  end
elseif ~(length(cfg.markersize)==1) && ~(length(cfg.markersize)==nchan)
  error('cfg.markersize should either have length = 1 or length = nchan')
end


% switch for brain plotting
if strcmp(cfg.renderbrain,'yes')
  % load template
  load(cfg.template)
  % exception
  if strcmp(cfg.template,'surface_white_both')
    cortex = [];
    cortex.tri  = bnd.tri;
    cortex.vert = bnd.pnt;
  end
  % Patch-reduction
  if ~strcmp(cfg.patchreduce,'no')
    % do it for the surface
    tempbrain.cortex.faces      = cortex.tri;
    tempbrain.cortex.vertices   = cortex.vert;
    tempbrain.cortex = reducepatch(tempbrain.cortex,cfg.patchreduce);
    cortex.tri  = tempbrain.cortex.faces;
    cortex.vert = tempbrain.cortex.vertices;
  end
  % plot brain itself
  hbrain = patch('faces',cortex.tri,'vertices',cortex.vert,'edgecolor','none','facecolor',cfg.braincolor); % ,[(1/255)*237 (1/255)*201 (1/255)*175]  [(1/255)*210 (1/255)*180 (1/255)*140]
  view(cfg.viewpoint)
  set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio', [1 1 1]);
  set(hbrain,'edgeAlpha',0)
  camlight(0,0);
  alpha(cfg.alpha);
  lighting phong
  material dull
  hold on
else
  view(cfg.viewpoint)
  set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio', [1 1 1]);
  hold on
end



% plot markers
camppos = get(gca,'cameraPosition');
if strcmp(cfg.plotelec,'yes')
  for ichan = 1:nchan
    %h = circlePlane3D(coords(ichan,:),camppos, cfg.markersize{ichan}, 0.1, 0, cfg.markercolor{ichan},'-k'); % circePlane3D DOESNT SAVE WELL TO EPS USING SAVEFIGSAMESIZE (not a circle, but connected points which are dpi dependent
    %set(h,'faceLighting','none')
    %set(h,'lineStyle','none')
    plot3(coords(ichan,1),coords(ichan,2),coords(ichan,3),cfg.markersymbol{ichan},'markerfacecolor',cfg.markercolor{ichan},'markeredgecolor',[0 0 0],'markersize',cfg.markersize{ichan});
  end
end
axis off
if ~strcmp(cfg.renderbrain,'yes')
  axis tight
end

% plot labels
if ~strcmp(cfg.label,'no')
  for ichan = 1:nchan
    text(coords(ichan,1),coords(ichan,2),coords(ichan,3),cfg.label{ichan},'interpreter','none','fontsize',cfg.labelsize,'color',cfg.labelcolor,'fontweight',cfg.labelweigth);
  end
end

% plot colorbar
if strcmp(cfg.colorbar,'yes')
  h = colorbar('location','NorthOutside');
  caxis([cmin cmax])
  colormap(cfg.colormap)
  pos = get(h,'position');
  pos(4) = pos(4)/3;
  pos(3) = (pos(3)/4)*3;
  axpos = get(gca,'outerposition');
  pos(2) = axpos(2);
  set(h,'position',pos)
end

% set axis limits
set(gca,'ylim',[-100 100])
set(gca,'xlim',[-100 100])
set(gca,'zlim',[-100 100])
if strcmp(cfg.axissqueeze,'yes')
  axis tight
end



% % delete
% if strcmp(cfg.renderbrain,'no')
%   set(gca,'xlimMode','manual')
%   set(gca,'ylimMode','manual')
%   set(gca,'zlimMode','manual')
%   delete(hbrain)
% end






























