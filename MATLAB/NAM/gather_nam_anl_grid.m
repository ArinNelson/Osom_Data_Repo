%=========================================================================%
% gather_nam_anl_grid.m
% Gather North American Mesoscale (NAM) model Analysis (ANL) grid data
% by Arin Nelson
% on 08/14/2021
%
% last edited by Arin Nelson on 08/14/2021
%=========================================================================%
addpath('../Utilities');
addpath('C:/Library/Codes/MATLAB/ROMS_matlab/m_map');

% Options
max_wait      = 0;	% max time to wait (secs) for web response (uses Java & parallelism, set to 0 to disable)
data_dir      = 'F:/OSOM_Data_Repo/NAM_ANL/'; % Location where to save
grid_filename = 'namanl_grid.nc';             % grid file name
lon_lim       = [-74.0 -69.0];                % longitude domain
lat_lim       = [ 40.0  42.5];                % latitude domain

% Constants
mask_url = 'https://www.ncei.noaa.gov/thredds/dodsC/model-namanl-old/201007/20100701/namanl_218_20100701_0000_000.grb';

%-------------------------------------------------------------------------%

% Full grid file directory
grid_file = [data_dir grid_filename];

% if the grid file already exists, ask user if they want to ovewrite it
check = false;
if(exist(grid_file,'file')==2)
	usr_input = input('A grid file with this name at this directory already exists.  Overwrite it? 1=yes, 0=no: ');
    if(usr_input==1)
        delete(grid_file);
    else
        check = true;
    end
end
if(check);  error('User aborted.'); end

% Read in additional grid variables
nam_x = ncread_web(max_wait,mask_url,'x');
nam_y = ncread_web(max_wait,mask_url,'y');

% x-y to lon-lat using m_map
clon      = -95.0;   %265.0;   
clat      = 25.0;
earth_rad = 6367.470;
[y,x]     = meshgrid(nam_y,nam_x);
m_proj('lambert conformal conic','clongitude',clon,'lat',[clat clat]);
[nam_lon,nam_lat] = m_xy2ll(x/earth_rad,y/earth_rad);

% find indices of grid over user-specified region
[ii,jj] = find(nam_lon>=lon_lim(1) & nam_lon<=lon_lim(end) & nam_lat>lat_lim(1) & nam_lat<lat_lim(end));
ii      = min(ii)-1 : max(ii)+1;
jj      = min(jj)-1 : max(jj)+1;

% These variables are what will be used by ncread, etc.
nam_i0 = ii(1);  nam_ni = numel(ii);
nam_j0 = jj(1);  nam_nj = numel(jj);
    
% Truncate grid
nam_lon  = nam_lon(ii,jj);
nam_lat  = nam_lat(ii,jj);
nam_x    = nam_x(ii);
nam_y    = nam_y(jj);
nam_mask = ncread_web(max_wait,mask_url,'Land_cover_land1_sea0_surface',[nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
nam_mask = 1-nam_mask;    % Want land=0, sea=1

% View
% imagesc(nam_x,nam_y,nam_mask'); set(gca,'ydir','normal');

%-------------------------------------------------------------------------%

% Generate grid file
nc_gen_nam_grid(grid_file,nam_ni,nam_nj);
  
% Save these variables to it
ncwrite(grid_file,'x',   nam_x   );
ncwrite(grid_file,'y',   nam_y   );
ncwrite(grid_file,'lon', nam_lon );
ncwrite(grid_file,'lat', nam_lat );
ncwrite(grid_file,'mask',nam_mask);
ncwrite(grid_file,'i0',  nam_i0  );
ncwrite(grid_file,'j0',  nam_j0  );

%=========================================================================%