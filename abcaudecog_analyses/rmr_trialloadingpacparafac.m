function nwaycomp = rmr_trialloadingpacparafac(nwaycomp,data,wplfcfg,compload)

% RMR_TRIALLOADINGPACPARAFAC
%
%  Use as:
%   nwaycomp = rmr_trialloadingpacparafac(cfg,nwaycomp,data)
%
% Obtain trial loadings for PARAFAC components extracted from 3/4-way wPLF arrays.
%
%
% Input:
%         wplfcfg = identical cfg as used rmr_phaseamplitudecoupling to obtain wPLFs that were subsequently decomposed
%        nwaycomp = output of PARAFAC decomposition on output of rmr_phaseamplitudecoupling without trial dimension
%            data = data structure as used for rmr_phaseamplitudecoupling
%
%
% Output:
%   nwaycomp = a data structure as obtained from nd_nwaydecomposition, as if a wPLF array having a trial dimension was used
%
%

% check 
if nargin==3
  compload = false;
else
  compload = istrue(compload);
end

% get ns
ntrial = numel(data.trial);
nchan  = numel(data.label);
hashfa = (isfield(wplfcfg,'amptstohfats') && istrue(wplfcfg.amptstohfats));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get trial loadings using NWAYDECOMP_PARAFAC snippets

% set size and number of modes and compflg (and default compmodes)
if hashfa
  nmode     = 4;
  smode     = [nchan nchan numel(wplfcfg.phasfreq) ntrial];
  compmodes = [1 1 0 compload];
else
  nmode = 5;
  smode = [nchan nchan numel(wplfcfg.ampfreq) numel(wplfcfg.phasfreq) ntrial];
  compmodes = [1 1 0 0 compload];
end
modeind = 1:nmode; % for small speedup when array is small

% get comp in parafac style
oldcomp = nwaycomp.comp;
ncomp   = numel(oldcomp);
comp = cell(1,nmode);
for imode = 1:nmode-1
  comp{imode} = zeros(smode(imode),ncomp);
  for icomp = 1:ncomp
    comp{imode}(:,icomp) = oldcomp{icomp}{imode};
  end
end
comp{end} = zeros(smode(end),ncomp);

% put back scaling
scaling = nwaycomp.scaling;
% magnitude scaling
comp{1} = comp{1} .* repmat(scaling{1},[smode(1) 1]);
% phase scaling
comp{1} = comp{1} ./ repmat(scaling{2},[smode(1) 1]);

% get trial loadings
for itrial = 1:ntrial
  disp(['--- obtaining trial loading for trial ' num2str(itrial) ' of ' num2str(ntrial) ' ---'])
  
  % set trial mode to current mode
  imode = modeind(end);
  
  % obtain current dat
  cfg = wplfcfg;
  cfg.trials = itrial;
  wplfdata = rmr_phaseamplitudecoupling(cfg,data);
  dat = squeeze(wplfdata.wplf);
  % unfold dat
  permorder   = [imode setdiff(modeind,imode)];
  reshapesize = [1 prod(smode(setdiff(modeind,imode)))];
  dat         = reshape(permute(dat,permorder),reshapesize);
 
  % set remaining mode, those to calculate the krb's over
  remmodes = modeind(modeind~=imode); % remaining modes
  
  % Calculate Z by a series of nested khatri-rao-bro products (using kr)
  Z = kr(comp(remmodes(end:-1:1)));

  % calculate ZctZ directly, which is faster than calculating Z first, especially when decomposing large arrays
  ZctZ = comp{remmodes(end)}'*comp{remmodes(end)};
  for iremmode = 1:numel(remmodes)-1
    ZctZ = ZctZ .* (comp{remmodes(iremmode)}'*comp{remmodes(iremmode)});
  end
  
  % Update the component matrix for the current mode
  % if complex, and currmode shouldn' be, take real of output that would otherwise be computed using catted real and imag parts, it's equivalent
  if compmodes(imode)==0 % computation below is identical (taking real of real data), but split up for code transparency
    %comp{imode} = real(dat{imode} * conj(Z)) * inv(real(ZctZ)).'; % real(Z'*X) = [Zre Zim]' * [Xre Xim];b  (.' is equal to conj in this case)
    comp{imode}(itrial,:) = real(dat * conj(Z)) / real(ZctZ);
  else
    %comp{imode} = (dat{imode} * conj(Z)) * conj(inv(ZctZ));
    comp{imode}(itrial,:) = (dat * conj(Z)) / conj(ZctZ);
  end

end


%%%% Post processing
% Magnitude postprocessing of components (make frobenius norm per loading vector norm = 1)
% make sure all explained variance is contained in first mode
for icomp = 1:ncomp
  % set mode1, which get's all the variance
  mode1 = comp{1}(:,icomp);
  % loop over the other modes
  for imode = 2:nmode
    currmode     = comp{imode}(:,icomp); % current mode of current compoment
    currmodenorm = norm(currmode,'fro');
    mode1        = mode1    .* currmodenorm;
    currmode     = currmode ./ currmodenorm;
    comp{imode}(:,icomp) = currmode;
  end
  % set mode1 back into original format
  comp{1}(:,icomp) = mode1;
end
disp('components have been magnitude normalized in all but the first mode')


% Phase postprocessing of complex components (make average magnitude weigthed phase per complex loading vector 0)
% make sure the scaling of phases is contained in first complex mode
compmodeindex = find(compmodes==1);
if sum(compmodes)>=2 % only do this if there is more than one complex mode
  for icomp = 1:ncomp
    % set mode1 (the first complex mode), which get's the phase scaling
    mode1 = comp{compmodeindex(1)}(:,icomp);
    % loop over the other modes
    for imode = compmodeindex(2:end)
      currmode   = comp{imode}(:,icomp); % current mode of current compoment
      meanangle  = angle(mean(currmode)); % angle of the mean weighted by magnitude
      phaseshift = exp(-1i*meanangle);
      mode1      = mode1    ./ phaseshift;
      currmode   = currmode .* phaseshift;
      comp{imode}(:,icomp) = currmode;
    end
    % set mode1 back into original format
    comp{compmodeindex(1)}(:,icomp) = mode1;
  end
  disp('components have been phase shifted so average magnitude-weighted-phase = 0 with respect to the first complex mode')
end
%%%% Post processing


%%%% Final Post processing
scaling = [];
% set norm of mode 1 per loading vector to 1 and save scaling coeffient
for icomp = 1:ncomp
  % set mode1
  mode1     = comp{1}(:,icomp);
  mode1norm = norm(mode1,'fro');
  mode1     = mode1 ./ mode1norm;
  comp{1}(:,icomp)  = mode1;
  scaling{1}(icomp) = mode1norm;
end
disp('first mode magnitude scaling coefficient removed')
% set average magnitude weighted phase to 0 for the first complex loading vector per component and save scaling coefficient
compmodeindex = find(compmodes==1);
if sum(compmodes)~=0 % only do this if there is a complex mode
  for icomp = 1:ncomp
    % set mode1
    mode1      = comp{compmodeindex(1)}(:,icomp); % compmodeindex determined above
    meanangle  = angle(mean(mode1)); % angle of the mean weighted by magnitude
    phaseshift = exp(-1i*meanangle);
    mode1      = mode1 .* phaseshift;
    comp{compmodeindex(1)}(:,icomp) = mode1;
    scaling{2}(icomp) = phaseshift;
  end
  disp('first complex mode has been phase shifted so average magnitude-weighted-phase = 0')
end
disp('post-processing finished')


% Calculate Tucker's Congruency coefficient
tuckcongr = ones(ncomp,ncomp);
for imode = 1:nmode
  innprod = comp{imode}' * comp{imode};
  tuckcongr = tuckcongr .* innprod;
end
% remove upper triangle and diagonal, make vector, and take abs
tuckcongr = tril(tuckcongr,-1);
tuckcongr(tuckcongr==0) = [];
tuckcongr = abs(tuckcongr);
if ncomp==1 % if ncomp is one, tuckers congruence cannot be calculated
  tuckcongr = NaN;
end
% display warning if too high
if max(tuckcongr) >= 0.85
  disp('Warning: some components are highly correlated, model might be degenerate')
end
%%%% Final Post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% reformat comp to nwaycomp format
outputcomp = [];
for icomp = 1:ncomp
  for idim = 1:numel(comp)
    outputcomp{icomp}{idim} = comp{idim}(:,icomp);
  end
end


% save in output nwaycomp
nwaycompout = nwaycomp;
nwaycompout.comp      = outputcomp;
nwaycompout.scaling   = scaling;
nwaycompout.tuckcongr = tuckcongr;
nwaycompout.dimord    = [nwaycompout.dimord '_rpt'];
nwaycomp = nwaycompout;




























