function [y,varargout] = roms_detide(x,t,tide_file)
%=========================================================================%
% y = roms_detide(x,t,tide_file)
% detide ROMS data x(t) using information on tides in tide_file
% 
% by Arin Nelson
% on 08/12/2021
%=========================================================================%
    
% FOR DEBUGGING
%tide_file = 'C:\Library\ROMS_Stuff\Resources\forcefiles_riroms\riroms_tides_2006.nc';

% Ensure missing values are NaN
x(abs(x)>1e30) = NaN;

% # dims in x
sizex = size(x);
n_dim = numel(sizex);

% Load tidal frequencies
if(exist(tide_file,'file')~=2)
    error(['File ' tide_file ' does not exist!']);
else
    tide_period = ncread(tide_file,'tide_period');
end

% Initialize trig args for tidal fits
tide_freq = 2*pi./(tide_period/24);  % Units of rads/day
n_tide    = numel(tide_freq);
n_time    = numel(t);
X         = NaN(n_time,2*n_tide+1);
X(:,1)    = 1;                      % For fitting the mean, should be 0 but just in case...
for iw=1:n_tide
      
	% Trig stuff  
    trigarg = tide_freq(iw).*t;
    cosvec  = cos( trigarg );
    sinvec  = sin( trigarg );
    
    % The X-values
    X(:,iw+1)        = cosvec;
    X(:,iw+1+n_tide) = sinvec;
    
    % Clean-up
    clear trigarg cosvec sinvec;
    
end
clear iw trigarg cosvec sinvec;

% Depending on if variable is 1D, 2D, or 3D, loop
y = NaN(size(x));
switch n_dim
    case 2
        T = NaN(size(x,1),n_tide);
        for ix=1:size(x,1)
            
            % This loc
        	xx = x(ix,:);   
            xx = xx(:);
            ii = find(~isnan(xx));
            if(~isempty(ii))
                
                % Tidal fit
                xx      = nandetrend(xx);
                [Fit,~] = lscov(X(ii,:),xx(ii),ones(numel(ii),1));
                
                % Create tidal time series
                tide_re  = Fit(          2:(n_tide+1) );
                tide_im  = Fit( (n_tide+2):end        );
                tide_fit = repmat(tide_re(:)' + sqrt(-1).*tide_im(:)',[nt 1]);
                tide_ts  = abs(tide_fit).*exp(sqrt(-1).*angle(tide_fit).*repmat(t(:),[1 n_tide]));
                
                % Subtract tidal time series
                y(ix,:) = x(ix,:)-sum(tide_ts,2);
                
                % Save tidal fit
                if(nargout>1); T(ix,:) = tide_re(:) + sqrt(-1).*tide_im(:); end
                
            end
            clear xx
            
        end
        clear ix;
    case 3
        T = NaN(size(x,1),size(x,2),n_tide);
        for ix=1:size(x,1)
        for iy=1:size(x,2)
            
            % This loc
        	xx = x(ix,iy,:);   
            xx = xx(:);
            ii = find(~isnan(xx));
            if(~isempty(ii))
                
                % Tidal fit
                xx      = nandetrend(xx);
                [Fit,~] = lscov(X(ii,:),xx(ii),ones(numel(ii),1));
                
                % Create tidal time series
                tide_re  = Fit(          2:(n_tide+1) );
                tide_im  = Fit( (n_tide+2):end        );
                tide_fit = repmat(tide_re(:)' + sqrt(-1).*tide_im(:)',[nt 1]);
                tide_ts  = abs(tide_fit).*exp(sqrt(-1).*angle(tide_fit).*repmat(t(:),[1 n_tide]));
                
                % Subtract tidal time series
                y(ix,iy,:) = x(ix,:)-sum(tide_ts,2);
                
                % Save tidal fit
                if(nargout>1); T(ix,iy,:) = tide_re(:) + sqrt(-1).*tide_im(:); end
                
            end
            clear xx
            
        end
        end
        clear ix iy;
    case 4
        T = NaN(size(x,1),size(x,2),size(x,3),n_tide);
        for ix=1:size(x,1)
        for iy=1:size(x,2)
        for iz=1:size(x,3)
            
            % This loc
        	xx = x(ix,iy,iz,:);   
            xx = xx(:);
            ii = find(~isnan(xx));
            if(~isempty(ii))
                
                % Tidal fit
                xx      = nandetrend(xx);
                [Fit,~] = lscov(X(ii,:),xx(ii),ones(numel(ii),1));
                
                % Create tidal time series
                tide_re  = Fit(          2:(n_tide+1) );
                tide_im  = Fit( (n_tide+2):end        );
                tide_fit = repmat(tide_re(:)' + sqrt(-1).*tide_im(:)',[nt 1]);
                tide_ts  = abs(tide_fit).*exp(sqrt(-1).*angle(tide_fit).*repmat(t(:),[1 n_tide]));
                
                % Subtract tidal time series
                y(ix,iy,iz,:) = x(ix,:)-sum(tide_ts,2);
                
                % Save tidal fit
                if(nargout>1); T(ix,iy,iz,:) = tide_re(:) + sqrt(-1).*tide_im(:); end
                
            end
            clear xx
            
        end
        end
        end
        clear ix iy iz;
    otherwise
        %???
end

end