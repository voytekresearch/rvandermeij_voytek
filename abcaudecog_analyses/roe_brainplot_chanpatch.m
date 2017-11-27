function roe_brainplot_chanpatch(cfg)

% ROE_BRAINPLOT_CHANPATCH plots channels on a talairach template brain
% It requires talairach coordinates, a viewpoint, and a set of channel loadings.
% Each 'patch' is the projected chanload value multiplied with a gaussian
% This function use the current figure window(!).
%
% The gaussian plotting code is taken from the LOC package made bij KJ Miller et al,
% described in Journal of Neuroscience Methods (2007).
%
% Use as
%    roe_brainplot_chanpatch(cfg)
%
%
% Necessary:
%    cfg.talairach     = chanX3 matrix of X,Y,Z coordinates
%    cfg.chanload      = chanX1 vector (if it is complex, only the angles will be used)
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
% 
% 
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
%    cfg.markerdisp    = marker-displacement in coordinates 
%    cfg.cparam        = vector to use for coloring markers (overwrites cfg.markercolor), should have dimensions channels*1
%    cfg.clim          = can be [cmin cmax], 'absmax' or 'minmax', must to be used when using cfg.cparam
%                        (default is 'absmax')
%    cfg.colormap      = colormap to be used when plotting patches and markers (when using 'cparam')
%    cfg.patchclim     = can be [cmin cmax], 'absmax' or 'minmax' (default is 'absmax')
%    cfg.colorbar      = 'yes' or 'no', displays colorbar
%    cfg.label         = 1xNchan cell-array or 'no'
%    cfg.labelsize     = fontsize
%    cfg.labelweigth   = font weight
%    cfg.labelcolor    = font color
%    cfg.patchreduce   = 'no', or a reduction-factor (e.g. 0.01)
%
%




% Set defaults
if ~isfield(cfg, 'alpha'),           cfg.alpha         = 1;                        end
if ~isfield(cfg, 'braincolor'),      cfg.braincolor    = [.9 .9 .9];               end
if ~isfield(cfg, 'markersize'),      cfg.markersize    = 3;                        end
if ~isfield(cfg, 'markercolor'),     cfg.markercolor   = 'g';                      end
if ~isfield(cfg, 'markersymbol'),    cfg.markersymbol  = '.';                      end
if ~isfield(cfg, 'markerdisp'),      cfg.markerdisp    = 0;                        end
if ~isfield(cfg, 'clim'),            cfg.clim          = 'absmax';                 end
if ~isfield(cfg, 'patchclim'),       cfg.patchclim     = 'absmax';                 end
if ~isfield(cfg, 'colormap'),        cfg.colormap      = 'jet';                    end
if ~isfield(cfg, 'colorbar'),        cfg.colorbar      = 'no';                     end
if ~isfield(cfg, 'label'),           cfg.label         = 'no';                     end
if ~isfield(cfg, 'labelsize'),       cfg.labelsize     = 8;                        end
if ~isfield(cfg, 'labelweigth'),     cfg.labelweigth   = 'normal';                 end
if ~isfield(cfg, 'labelcolor'),      cfg.labelcolor    =  [1 1 1];                 end
if ~isfield(cfg, 'patchreduce'),     cfg.patchreduce   = 'no';                     end
if ~isfield(cfg, 'viewpoint'),       cfg.viewpoint     = [90 0];                   end
if ~isfield(cfg, 'displace'),        cfg.displace      = [0 0 0];                  end
if ~isfield(cfg, 'talairach'),       error('talairach coordinates needed'),        end



% set nchan and talcoords
nchan = size(cfg.talairach,1);
talcoords = cfg.talairach;


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
  colormat = eval([cfg.colormap '(10e6);']);
  
  % Transform cparam to go from 0 1 and multiply with 10e6 to act as index into colormat
  cparam = round(((cfg.cparam + -cmin) .* (1 / (-cmin + cmax))) .* 10e6);
  cparam(cparam==0) = 1; % if a color index is zero, zet it to the first color-value (which has index 1)
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


%load('tal_whole.mat')
load('wholebrain.mat')
% Patch-reduction
if ~strcmp(cfg.patchreduce,'no')
  % do it for the surface
  tempbrain.cortex.faces      = cortex.tri;
  tempbrain.cortex.vertices   = cortex.vert;
  tempbrain.cortex = reducepatch(tempbrain.cortex,cfg.patchreduce);
  cortex.tri  = tempbrain.cortex.faces;
  cortex.vert = tempbrain.cortex.vertices;
end



% prepare patch
chanload = cfg.chanload;
% build distance matrix to use for interpolation limits
distmat = NaN(nchan,nchan);
for iseed = 1:nchan
  xs = talcoords(iseed,1);
  ys = talcoords(iseed,2);
  zs = talcoords(iseed,3);
  for itarg = 1:nchan
    xt = talcoords(itarg,1);
    yt = talcoords(itarg,2);
    zt = talcoords(itarg,3);
    % Calculate Euclidian distance
    distmat(iseed,itarg) = ((xs-xt)^2 + (ys-yt)^2 + (zs-zt) ^2 ) ^ 0.5;
  end
end
distmat(logical(eye(nchan))) = inf;
patchdistlim = mean(min(distmat));
% if complex, to be plotted colors are angles (so normalize to unit)
if ~isreal(chanload)
  chanload = chanload ./ abs(chanload);
end

% build face-vertex-color data
col = zeros(size(cortex.vert,1),1);
gausspread = 50; % in mm
for ichan = 1:length(talcoords(:,1))
  bx = abs(cortex.vert(:,1) - talcoords(ichan,1));
  by = abs(cortex.vert(:,2) - talcoords(ichan,2));
  bz = abs(cortex.vert(:,3) - talcoords(ichan,3));
  %d = chanload(ichan) .* (1/(bx + bz + by))
  d = chanload(ichan) * exp((-(bx .^ 2 + bz .^ 2 + by .^ 2)) / gausspread);
  col = col + (d ./ nchan);
end
% if complex, to be plotted colors are angles
if ~isreal(col)
  col = angle(col);
end

% index colors and apply patchdistlim
facecol = repmat(cfg.braincolor,size(cortex.vert,1),1);
colormat = eval([cfg.colormap '(10e6);']);
if strcmp(cfg.patchclim,'maxmin')
  cmin = min(col);
  cmax = max(col);
elseif strcmp(cfg.patchclim,'absmax')
  cmin = -max(max(abs(col)));
  cmax = max(max(abs(col)));
else
  cmin = cfg.patchclim(1);
  cmax = cfg.patchclim(2);
end
for ivert = 1:size(facecol,1)
  % detect whether it falls within patchdistlim
  if all(sum(abs(repmat(cortex.vert(ivert,:),[nchan,1]) - talcoords) < patchdistlim))
    % index the color from col
    currcol = col(ivert);
    currcol = round(((currcol + -cmin) .* (1 / (-cmin + cmax))) .* 10e6);
    if currcol==0
      currcol = 1; % if color index is zero, set it to the first color-value (which has index 1)
    end
    facecol(ivert,:) = colormat(currcol,:);
  end
end

% plot brain surface
trisurf(cortex.tri, cortex.vert(:, 1), cortex.vert(:, 2), cortex.vert(:, 3), 'FaceVertexCData', facecol);
view(cfg.viewpoint)
hold on
axis off
set(gcf,'Renderer', 'zbuffer')
set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio', [1 1 1]);
lighting gouraud;
material dull;
alpha(cfg.alpha);
camlight(0,0);
shading interp;
colormap(hsv)


% displace talcoords before plotting markers
talcoords = talcoords + repmat(cfg.displace,[size(talcoords,1) 1]);

% plot markers
for ichan = 1:nchan
  plot3(talcoords(ichan,1),talcoords(ichan,2),talcoords(ichan,3),cfg.markersymbol{ichan},'Color',cfg.markercolor{ichan},'markersize',cfg.markersize{ichan},'LineWidth',3);
end
axis off

% plot labels
if ~strcmp(cfg.label,'no')
  for ichan = 1:nchan
    text(talcoords(ichan,1),talcoords(ichan,2),talcoords(ichan,3),cfg.label{ichan},'interpreter','none','fontsize',cfg.labelsize,'color',cfg.labelcolor,'fontweight',cfg.labelweigth);
  end
end

% plot colorbar
if strcmp(cfg.colorbar,'yes')
  h = colorbar('location','NorthOutside');
  caxis([cmin cmax])
  colormap(cfg.colormap)
  pos = get(h,'position');
  pos(4) = pos(4)/4;
  pos(3) = (pos(3)/4)*3;
  set(h,'position',pos)
end
  

