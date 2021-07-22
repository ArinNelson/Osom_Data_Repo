%=========================================================================%
% gather_nam_anl_grid.m
% Gather North American Mesoscale (NAM) model Analysis (ANL) grid data
% by Arin Nelson
% on 07/21/2021
%
% Example file url for testing/debugging:
% dat_url = 'https://www.ncei.noaa.gov/thredds/dodsC/model-namanl-old/201001/20100101/namanl_218_20100101_0000_000.grb';
%
% last edited by Arin Nelson on 07/21/2021
%=========================================================================%

% Options
max_wait        = 120;	% max time to wait (secs) for web response (uses Java & parallelism, set to <=0 to disable)
grid_savedir    = 'F:/OSOM_Data_Repo/NAM/';     % location to save grid file
grid_filename   = [data_dir 'nam_grid.nc'];     % grid file name
lon_lim         = [-74.0 -69.0];                % longitude domain
lat_lim         = [ 40.0  42.5];                % latitude domain

% Constants
grid_url = 'https://rda.ucar.edu/datasets/ds609.0/docs/latlon-g218.txt';

%-------------------------------------------------------------------------%

% Full grid file directory
grid_file = [grid_savedir grid_filename];

% if the grid file already exists, ask user if they want to ovewrite it
check = false;
if(exist(grid_file,'file')==2
	usr_input = input('A grid file with this name at this directory already exists.  Overwrite it? 1=yes, 0=no');
    if(usr_input==1)
        delete(grid_file);
    else
        check = true;
    end
end
if(check);  error('User aborted.'); end

% read in NAM 218 lon/lat grid (from url)
txt_raw = textscan( webread(grid_url,weboptions('ContentType','text')) , '%s', 'delimiter', '\n' );
txt_raw = txt_raw{1};
txt_raw(1:2) = [];
txt_parse = zeros(4,numel(txt_raw));
for i=1:numel(txt_raw)
	tmp = textscan(txt_raw{i},'%f %f %f %f');
	txt_parse(:,i) = [tmp{:}];
end
clear i tmp txt_raw;

% convert to 2D grid
nx      = max(txt_parse(1,:));
ny      = max(txt_parse(2,:));
nam_lon = zeros(nx,ny);
nam_lat = zeros(nx,ny);
for i=1:size(txt_parse,2)
	nam_lat(txt_parse(1,i),txt_parse(2,i)) =  txt_parse(3,i);
	nam_lon(txt_parse(1,i),txt_parse(2,i)) = -txt_parse(4,i); % Units are in deg.W, but want deg.E, so set to negative
end
clear txt_parse i nx ny;

% find indices of grid over user-specified region
[ii,jj] = find(nam_lon>=lon_lim(1) & nam_lon<=lon_lim(end) & nam_lat>lat_lim(1) & nam_lat<lat_lim(end));
ii      = min(ii)-1 : max(ii)+1;
jj      = min(jj)-1 : max(jj)+1;

% These variables are what will be used by ncread, etc.
nam_i0 = ii(1);  nam_ni = numel(ii);
nam_j0 = jj(1);  nam_nj = numel(jj);
    
% Truncate grid
nam_lon = nam_lon(ii,jj);
nam_lat = nam_lat(ii,jj);

% Read in additional grid variables
nam_x    = ncread_web(max_wait,mask_url,'x',nam_j0,nam_ni);
nam_y    = ncread_web(max_wait,mask_url,'y',nam_j0,nam_nj);
nam_mask = ncread_web(max_wait,mask_url,'Land_cover_land1_sea0_surface',[nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
nam_mask = 1-nam_mask;    % Want land=0, sea=1

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