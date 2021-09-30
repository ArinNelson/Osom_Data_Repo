function varargout = nanfill(varargin)
%--------------------------------------------------------------------------
% nanfill - 1-D interpolation of NaN values
% by Arin Nelson
% last updated 09/30/2021
% 
% xf = nanfill(x)
% fills in NaNs in x using linear interpolation and outputs result in xf
% 
% xf = nanfill(t,x)
% fills in NaNs in x(t) using linear interpolation and outputs result in xf
% note: t must be in ascending or descending order
% 
% [tf,xf] = nanfill(t,x)
% fills in NaNs in t using linear interpolation, then fills in x(t)
% 
% [...] = nanfill(...,'opt_name',opt_value)
% specifies any number of options:
% 
%   'method': specify string of method to use in interp1 function
%             default: 'linear'
%             type 'help interp1' to see a list of available methods
%   
%--------------------------------------------------------------------------

    % Set option defaults
    method = 'linear';

    % Parse input data
    if(nargin>0)
        if(mod(nargin,2)==1)
            x = varargin{1};
            t = 1:numel(x);
        elseif(mod(nargin,2)==0)
            t = varargin{1};
            x = varargin{2};
        end
    else
        error('nanfill requires at least 1 input argument!');
    end
    
    % Parse input options
    if(nargin>2)
        if(mod(nargin,2)==1)
            for i=2:2:nargin-1
                eval([varargin{i} '=' varargin{i+1} ';']);
            end
        elseif(mod(nargin,2)==0)
            for i=3:2:nargin-1
                eval([varargin{i} '=' varargin{i+1} ';']);
            end
        end
    end
    
    % Do interpolation of missing times
    ii = find(isnan(t));
    if(~isempty(ii))
        jj    = find(~ismember(1:numel(x),ii)==1);
        t(ii) = interp1(jj,t(jj),ii,method);
    end
        
    % Do interpolation of missing data points
    ii = find(isnan(x));
    if(~isempty(ii))
        jj    = find(~ismember(1:numel(x),ii)==1);
        x(ii) = interp1(t(jj),x(jj),t(ii),method);
    end
    
    % Return requested outputs
    switch nargout
        case 1;     varargout = {x};
        case 2;     varargout = {t,x};
    end
    
end