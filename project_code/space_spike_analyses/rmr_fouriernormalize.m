function [fourier, scaling] = rmr_fouriernormalize(fourier,normmethod)

% RMR_FOURIERNORMALIZE is a simple function used for normalizing the output of the rmr_SPACEinput_XXX functions.
% Not wonderfully documented, use at your own risk... (i.e. go carefully the operations you're requesting)
% 
%  Use as:
%  fourier = rmr_fouriernormalize(fourier,normmethod)
%
%  The following normalizations are available:
%               (see the code below for their details)
%          'avgoverepoch'
%          'coh'
%          'epochcoh'
%          'powerepochnorm'
%          'epochkthrootpower'
%          'offdiagpower'
%          'propodtod'
%          'kthrootpower'  (specify as '2throotpower','8throotpower', '2048throotpower', etc)
%          'reusescaling'  (specify normmethod as the 1xNchan scaling array received as output prevously)
%                          (only applicable for scaling involving a single coefficient per channel)
%          'none' 
%          
%  Normalizations can be performed in sequence using underscores. For example, 'avgoverepoch_8throotpower'
%  first performs 'avgoverepoch' and then on its result '8throotpower'.
%
%
%       Input:
%         fourier = 4-way chan_freq_epoch_taper array of fourier coefficients
%
%      Output:
%         fourier = 4-way chan_freq_epoch_taper array of fourier coefficients (SINGLE precision)
%
%
% Copyright (C) 2016-present, Roemer van der Meij


% parse input
if isnumeric(normmethod)
  scaling    = normmethod;
  normmethod = 'reusescaling';
end

% parse multiple options using '_'
normmethod = tokenize(normmethod,'_');
for itok = 1:(numel(normmethod)-1)
  fourier = roe_fouriernormalize(fourier,normmethod{itok});
end
normmethod = normmethod{end};
  
% parse kthrootpower
if numel(normmethod)>=12 && strcmp(normmethod(end-10:end),'throotpower') && ~strcmp(normmethod(1:5),'epoch')
  k          = str2double(normmethod(1:end-11));
  normmethod = 'kthrootpower';
end
% parse epochkthrootpower
if numel(normmethod)>=12 && strcmp(normmethod(end-10:end),'throotpower') && strcmp(normmethod(1:5),'epoch')
  k          = str2double(normmethod(6:end-11));
  normmethod = 'epochkthrootpower';
end


% switch over normalization methods
switch normmethod
  
  
  case 'avgoverepoch'
    %%%%%%%%%%%%
    % average over epochs
    % replace the epoch specific fourier coefficients, for each frequency,
    % by their average over epochs
    fourold = fourier;
    fourier = complex(NaN(size(fourold,1),size(fourold,2),1,size(fourold,1)),NaN(size(fourold,1),size(fourold,2),1,size(fourold,1)));
    for ifreq = 1:size(fourold,2)
      csd = zeros(size(fourold,1));
      for iepoch = 1:size(fourold,3)
        currfour  = double(squeeze(fourold(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        csd = csd + ((currfour*currfour') ./ size(fourold,3));
      end
      
      % Reduce SPACE memory load and computation time by replacing each chan_taper matrix by the
      % Eigenvectors of its chan_chan cross-products weighted by sqrt(Eigenvalues).
      % This possible because (1) SPACE only uses the cross-products of the chan_taper matrices
      % (i.e. the frequency- and epoch-specific CSD) and (2) the Eigendecomposition of a symmetric
      % matrix A is A = VLV'.
      % As such, VL^.5 has the same cross-products as the original chan_taper matrix.
      [V L] = eig(csd);
      L     = diag(L);
      tol   = max(size(csd))*eps(max(L)); % compute tol using matlabs default
      zeroL = L<tol;
      eigweigth = V(:,~zeroL)*diag(sqrt(L(~zeroL)));
      % positive semi-definite failsafe
      if any(L<-tol)
        error('csd not positive semidefinite')
      end
      % save in fourier
      currm = size(eigweigth,2);
      fourier(:,ifreq,1,1:currm) = eigweigth;
    end
    % trim fourier (ntaper is a bit unpredictable, hence the need for trimming)
    notnan = logical(squeeze(sum(~isnan(squeeze(fourier(1,:,:,:))),1)));
    fourier = fourier(:,:,:,notnan);
    %%%%%%%%%%%%
    
   

  case 'coh'
    %%%%%%%%%%%%
    % transform fourier into 'coherency'
    % i.e., scale fourier such that the power, of each channel,
    % averaged over frequences and epochs, is equal to 1
    power = zeros(size(fourier,1),size(fourier,2));
    for iepoch = 1:size(fourier,3)
      currfour = double(squeeze(fourier(:,:,iepoch,:)));
      power = power + nansum(abs(currfour).^2,3);
    end
    power = power ./ size(fourier,3);
    for iepoch = 1:size(fourier,3)
      for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        currfour  = bsxfun(@rdivide,currfour,sqrt(power(:,ifreq)));
        % save in fourier
        fourier(:,ifreq,iepoch,nonnanind) = single(currfour);
      end
    end
    %%%%%%%%%%%%
    
   
    
  case 'epochcoh'
    %%%%%%%%%%%%
    % transform fourier into 'coherency', for every epoch
    % i.e., scale fourier such that the power, of each channel
    % FOR EACH epoch/freq is 1.
    for iepoch = 1:size(fourier,3)
      for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        currpower = nansum(abs(currfour).^2,2);
        actind    = currpower~=0; % indices of channels to normalize
        currfour(actind,:) = bsxfun(@rdivide,currfour(actind,:),sqrt(currpower(actind)));
        
        % save in fourier
        fourier(:,ifreq,iepoch,nonnanind) = single(currfour);
      end
    end
    %%%%%%%%%%%%

   
    
    
  case 'powerepochnorm'
    %%%%%%%%%%%%
    % scale fourier such that the power, of each channel, 
    % FOR EACH epoch/freq of the CSD, is equal to that of its average
    power = zeros(size(fourier,1),size(fourier,2));
    for iepoch = 1:size(fourier,3)
      currfour = double(squeeze(fourier(:,:,iepoch,:)));
      power = power + nansum(abs(currfour).^2,3);
    end
    power = power ./ size(fourier,3);
    % scale power by kth root
    scaling = power;
    for iepoch = 1:size(fourier,3)
      for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        currpower = nansum(abs(currfour).^2,2);
        actind    = currpower~=0; % indices of channels to normalize
        currfour(actind,:) = bsxfun(@rdivide,currfour(actind,:),sqrt(currpower(actind)));
        currfour  = bsxfun(@times,currfour,sqrt(scaling(:,ifreq)));
        % save in fourier
        fourier(:,ifreq,iepoch,nonnanind) = single(currfour);
      end
    end
    %%%%%%%%%%%%
    
    
    
  case 'epochkthrootpower'
    %%%%%%%%%%%%
    % scale fourier such that the power, of each channel, 
    % FOR EACH epoch/freq of the CSD, is equal to that of its kth root
    for iepoch = 1:size(fourier,3)
      for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        currpower = nansum(abs(currfour).^2,2);
        actind    = currpower~=0; % indices of channels to normalize
        currfour(actind,:) = bsxfun(@rdivide,currfour(actind,:),sqrt(currpower(actind)));
        currfour  = bsxfun(@times,currfour,sqrt(currpower(actind).^(1/k)));
        % save in fourier
        fourier(:,ifreq,iepoch,nonnanind) = single(currfour);
      end
    end
    %%%%%%%%%%%%
    
    
    
  case 'offdiagpower'
    %%%%%%%%%%%%
    % scale fourier such that the average power, of each channel,
    % over epochs/freq of the CSD, is equal to the mean of the abs of the off-diagonal
    % elements (abs and mean over freq)
    csd = zeros(size(fourier,1),size(fourier,1),size(fourier,2));
    for iepoch = 1:size(fourier,3)
      for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        csd(:,:,ifreq) = csd(:,:,ifreq) + ((currfour * currfour') ./ size(fourier,3));
      end
    end
    csd   = mean(abs(csd),3);
    power = diag(csd);
    covar = sum(csd - diag(power),2) ./ (size(fourier,1)-1);
    % combine
    scaling = covar ./ power;
    for iepoch = 1:size(fourier,3)
      for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        currfour  = bsxfun(@times,currfour,sqrt(scaling));
        % save in fourier
        fourier(:,ifreq,iepoch,nonnanind) = single(currfour);
      end
    end
    %%%%%%%%%%%%
    
   
    
  case 'propodtod'
    %%%%%%%%%%%%
    % scale fourier such that the average power, of each channel,
    % over epochs/freq of the CSD, is equal to the proportion of its power to
    % the mean of the abs of the off-diagonal elements (abs and mean over freq)
    csd = zeros(size(fourier,1),size(fourier,1),size(fourier,2));
    for iepoch = 1:size(fourier,3)
      for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        csd(:,:,ifreq) = csd(:,:,ifreq) + ((currfour * currfour') ./ size(fourier,3));
      end
    end
    csd   = mean(abs(csd),3);
    power = diag(csd);
    covar = sum(csd - diag(power),2) ./ (size(fourier,1)-1);
    propodtod = covar ./ power;
    % combine
    scaling = propodtod ./ power;
    for iepoch = 1:size(fourier,3)
      for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        currfour  = bsxfun(@times,currfour,sqrt(scaling));
        % save in fourier
        fourier(:,ifreq,iepoch,nonnanind) = single(currfour);
      end
    end
    %%%%%%%%%%%%
    
    
  case 'kthrootpower'
    %%%%%%%%%%%%
    % scale fourier such that the average power, of each channel, 
    % over epochs/freq of the CSD, is equal to its kth root
    power = zeros(size(fourier,1),1);
    for iepoch = 1:size(fourier,3)
      currfour = double(squeeze(fourier(:,:,iepoch,:)));
      power = power + nansum(nansum(abs(currfour).^2,3),2);
    end
    power = power ./ (size(fourier,3) .* size(fourier,2));
    % scale power by its kth root
    scaling = power .^ (1/k);
    for iepoch = 1:size(fourier,3)
      for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        currfour  = bsxfun(@rdivide,currfour,sqrt(power));
        currfour  = bsxfun(@times,currfour,sqrt(scaling));
        % save in fourier
        fourier(:,ifreq,iepoch,nonnanind) = single(currfour);
      end
    end
    %%%%%%%%%%%%
 
    
  case 'reusescaling'
    %%%%%%%%%%%%
    % scale with previous scaling coefficients given as input
    % (only applicable to scaling involving a single coefficient per channel)
    power = zeros(size(fourier,1),1);
    for iepoch = 1:size(fourier,3)
      currfour = double(squeeze(fourier(:,:,iepoch,:)));
      power = power + nansum(nansum(abs(currfour).^2,3),2);
    end
    power = power ./ (size(fourier,3) .* size(fourier,2));
    for iepoch = 1:size(fourier,3)
      for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        currfour  = bsxfun(@rdivide,currfour,sqrt(power));
        currfour  = bsxfun(@times,currfour,sqrt(scaling));
        % save in fourier
        fourier(:,ifreq,iepoch,nonnanind) = single(currfour);
      end
    end
    %%%%%%%%%%%%
       
    
    
  case 'none'
    % do nothing
    
  otherwise
    error('normmethod not implemented')
end





 
  

    
    




