function [var_value, varargout] = roms_gather_output_var(roms_dir, var_name, varargin)
%=========================================================================%
% [var_value, time] roms_gather_output_var(roms_dir, var_name)
% Gather variable defined by string var_name from the ROMS output files
% prefixed by the string roms_dir.
% 
% Example:
% [zeta, time] = roms_gather_output_var('/ROMS/roms_his_','zeta');
%=========================================================================%
% by Arin Nelson 
% on 10/20/2020
% 
% last updated by Arin Nelson on 08/10/2021
%=========================================================================%

% Parse inputs
if(nargin>2)
    depth_level = varargin{1};
else
    depth_level = 1;
end

% Gather list of files
file_list = ls([roms_dir '*.nc']);

% Gather variable information from first data file
file_info = ncinfo([output_dir file_list(1,:)]);
var_info  = file_info.Variables( strcmp({file_info.Variables.Name},var_name)==1 );
if(isempty(var_info)); error(['Variable ' var_name ' not found in ROMS output!']); end

% Determine if variable is 2D or 3D
var_dims = numel(var_info.Size)-1;

% Get # times
n_files     = size(file_list,1);
nt_per_file = var_info.Size(end);
nt          = n_files*nt_per_file;

% Init time variable if wanted
if(nargout>1)
    time = zeros(nt,1);	
end

% Initialize data variable
var_value = NaN(var_info.Size(1),var_info.Size(2),nt);

% Loop through files
for i=1:n_files

	% Time indices
    i1 = nt_per_file*(i-1) + 1;
    i2 = nt_per_file*i;
    ii = i1:i2;
  
    % Load time
    if(nargout>1)
        time(ii) = ncread([roms_dir file_list(i,:)],'ocean_time');
    end
  
    % Load data
    if(var_dims==2)
        var_value(:,:,ii) = ncread([output_dir file_list(i,:)],var_name);
    else
        var_value(:,:,ii) = squeeze( ncread([output_dir file_list(i,:)],var_name,[1 1 depth_level 1],[inf inf 1 inf]) );
    end
    
    % Clean-up
    clear i1 i2 ii;
    
end

end