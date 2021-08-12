%=========================================================================%
% gather_narr.m
% Gather North American Regional Reanalysis (NARR) data and sort
% by Arin Nelson
% on 07/17/2021
%
% Catalog of data
% dat_url = 'https://psl.noaa.gov/thredds/catalog/Datasets/NARR/monolevel/catalog.html';
% 
% Example file for debugging
% https://psl.noaa.gov/thredds/dodsC/Datasets/NARR/monolevel/pottmp.hl1.1979.nc
%
% last edited by Arin Nelson on 07/17/2021
%=========================================================================%
clc; clear mex; addpath('../Utilities');

% Switches
Switch    = zeros(9,1);
Switch(1) = 0;      % Initialize grid file and sampling info
Switch(2) = 1;      % Perform data gathering

% Options
date_start = [2018,01];     % Start year, month to gather data for
date_end   = [2021,01];     % End   year, month to gather data for
max_wait   = 120;           % Max time to wait (secs) for web response
%data_dir   = '/gpfs/data/epscor/anelson5/OSOM_Data_Repo/NARR/';
data_dir   = 'F:/OSOM_Data_Repo/NARR/';
grid_file  = [data_dir 'narr_grid.nc'];
var_to_get = {'winds','Pair','Tair','Qair','rain','lwrad','swrad'};

%-------------------------------------------------------------------------%

% Constants
date_min = [1979,01,01];
date_max = [2020,12,31];
base_url = 'https://psl.noaa.gov/thredds/dodsC/Datasets/NARR/monolevel/';
grid_url = 'https://psl.noaa.gov/thredds/dodsC/Datasets/NARR/monolevel/pottmp.hl1.1979.nc'; 
lon_lim  = [-74.0 -69.0];
lat_lim  = [ 40.0  42.5];
n_var    = numel(var_to_get);

% Variable info
% List of variables available: https://psl.noaa.gov/data/gridded/data.narr.monolevel.html
var_info = {'Uwind',        'wind_time',    'uwnd.10m',             '';         ...	% units of m/s
            'Vwind',        'wind_time',    'vwnd.10m',             '';         ...	% units of m/s
            'Pair',         'pair_time',    'prmsl',                './100';    ...	% Convert Pa to mbar
            'Tair',         'tair_time',    'pottmp.sfc',           '-273.15';  ...	% Convert K to C
            'Qair',         'qair_time',    'rhum.2m',              '';         ...	% units of percent
            'Cfra',         'cloud_time',   'cdlyr',                './100';    ...	% units of percent to fraction
            'rain',         'rain_time',    'prate',             	'';         ...	% Units of kg/m2/s
            'lwrad_down',   'lrf_time',     'dlwrf',                '';         ... % Units of W/m2
            'lwrad_up',     'lrf_time',     'ulwrf',                '';         ... % Units of W/m2
            'swrad_down',   'srf_time',     'dswrf',                '';         ... % Units of W/m2
            'swrad_up',     'srf_time',     'uswrf',                '';         ... % Units of W/m2
            'Tdew',         'tdew_time',    'dpt.2m',               '-273.15';  ... % Convert K to C
            'sensible',     'shf_time',     'shtlf',                '';         ... % Units of W/m2
            'latent',       'lhf_time',     'lhtlf',                '';         ... % Units of W/m2
           };

%=========================================================================%
if(Switch(1)==1)
    
  % Delete grid file if it already exists  
  if(exist(grid_file,'file')==2)
    delete(grid_file);
  end
  
  % Read in NARR lon/lat grid (from url)
  narr_x   = ncread_web(max_wait,grid_url,'x'  );
  narr_y   = ncread_web(max_wait,grid_url,'y'  );
  narr_lon = ncread_web(max_wait,grid_url,'lon');
  narr_lat = ncread_web(max_wait,grid_url,'lat');
    
  % Find grid over specified region
  [ii,jj] = find(narr_lon>=lon_lim(1) & narr_lon<=lon_lim(end) & narr_lat>lat_lim(1) & narr_lat<lat_lim(end));
  ii      = min(ii)-1 : max(ii)+1;
  jj      = min(jj)-1 : max(jj)+1;
  
  % These variables are what will be used by ncread, etc.
  narr_i0 = ii(1);  ni = numel(ii);
  narr_j0 = jj(1);  nj = numel(jj);
    
  % Truncate grid
  narr_x   = narr_x(ii);
  narr_y   = narr_y(jj);
  narr_lon = narr_lon(ii,jj);
  narr_lat = narr_lat(ii,jj);
  
  % Generate grid file
  nc_gen_narr_grid(grid_file,ni,nj);
    
  % Save these variables to it
  ncwrite(grid_file,'x',   narr_x   );
  ncwrite(grid_file,'y',   narr_y   );
  ncwrite(grid_file,'lon', narr_lon );
  ncwrite(grid_file,'lat', narr_lat );
  ncwrite(grid_file,'i0',  narr_i0  );
  ncwrite(grid_file,'j0',  narr_j0  );
  
  % Land-sea mask (gathered manually from Tair for 2010/01/01 00:00:00, 
  % as there's a ~5C difference between land and water at this time
  load('narr_mask.mat');
  ncwrite(grid_file,'mask',int32(narr_mask));
  
end
%=========================================================================%
if(Switch(2)==1)

  % Load necessary variables
  narr_mask = ncread(grid_file,'mask');
  narr_i0   = ncread(grid_file,'i0');
  narr_j0   = ncread(grid_file,'j0');
  [narr_ni,narr_nj] = size(narr_mask);

  
  % Save directories
  for i=1:numel(var_to_get)
    this_dir = [data_dir var_to_get{i}];
    if(exist(this_dir,'dir')~=7)
      mkdir(this_dir);
    end
  end
  
  % Start year and month
  year_on  = date_start(1);
  month_on = date_start(2);
  
  % Loop through available months and years
  while( (year_on + (month_on-0.5)/12) < (date_end(1) + (date_end(2)-0.25)/12) )
    
    % Time index for beginning of month in input file
    if(month_on==1);    i0_thismonth = 1;
    else;               i0_thismonth = sum(eomday(year_on,1:month_on-1))*8+1;
    end;                n_thismonth  = sum(eomday(year_on,1:month_on))*8 - i0_thismonth+1;
    
    % Time vector for this output file (units of days since beginning of month)
    t = 0 : 0.125 : (n_thismonth/8-1)+0.8751;
    
    % Initialize save files for this year & month
    save_file = cell(n_var,1);
    for iv=1:n_var
      save_file{iv} = [data_dir var_to_get{iv} '/NARR_' var_to_get{iv} '_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
      if(exist(save_file{iv},'file')~=2)
        nc_gen_narr_data(save_file{iv},narr_ni,narr_nj,n_thismonth,var_to_get(iv));  % NOTE PARENTHESIS ON var_to_get(iv)!
      end
    end
    clear i
    
    % Loop through variables
    clc; disp(['Gathering NARR data for year/month ' num2str(year_on) '/' sprintf('%0.2d',month_on) '...']);
    for iv=1:n_var
    switch var_to_get{iv}
       
      %-----------------------------------------------------------------%
      case 'winds'
            
         u_url  = [base_url 'uwnd.10m.' num2str(year_on) '.nc']; 
         v_url  = [base_url 'vwnd.10m.' num2str(year_on) '.nc']; 
         narr_u = ncread(u_url,'uwnd',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         narr_v = ncread(v_url,'vwnd',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         ncwrite(save_file{iv},'Uwind',narr_u);
         ncwrite(save_file{iv},'Vwind',narr_v);
         ncwrite(save_file{iv},'wind_time',t);
         clear u_url v_url narr_u narr_v;
            
      %-----------------------------------------------------------------%    
      case 'Pair'
        
         z_url  = [base_url 'prmsl.' num2str(year_on) '.nc'];
         narr_z = ncread(z_url,'prmsl',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]);
         narr_z = narr_z ./ 100;
         ncwrite(save_file{iv},'Pair',narr_z);
         ncwrite(save_file{iv},'pair_time',t);
         clear z_url narr_z;
         
      %-----------------------------------------------------------------%    
      case 'Tair'
        
         z_url  = [base_url 'pottmp.sfc.' num2str(year_on) '.nc'];
         narr_z = ncread(z_url,'pottmp',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         narr_z = narr_z - 273.15;
         ncwrite(save_file{iv},'Tair',narr_z);
         ncwrite(save_file{iv},'tair_time',t);
         clear z_url narr_z;
         
      %-----------------------------------------------------------------%    
      case 'Qair'
        
         z_url  = [base_url 'rhum.2m.' num2str(year_on) '.nc'];
         narr_z = ncread(z_url,'rhum',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         ncwrite(save_file{iv},'Qair',narr_z);
         ncwrite(save_file{iv},'qair_time',t);
         clear z_url narr_z;
         
      %-----------------------------------------------------------------%    
      case 'Cfra'
        
         %z1_url  = [base_url 'cdcon.' num2str(year_on) '.nc'];
         z2_url  = [base_url 'cdlyr.' num2str(year_on) '.nc'];
         %narr_z1 = ncread(z1_url,'cdcon',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         narr_z2 = ncread(z2_url,'cdlyr',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         narr_z  = narr_z2;         % + narr_z1;
         narr_z  = narr_z./100;     % Convert % to fraction
         ncwrite(save_file{iv},'Cfra',narr_z);
         ncwrite(save_file{iv},'cloud_time',t);
         clear z2_url narr_z2 narr_z; % z1_url narr_z1; 
         
      %-----------------------------------------------------------------%    
      case 'rain'
        
         z_url  = [base_url 'prate.' num2str(year_on) '.nc'];
         narr_z = ncread(z_url,'prate',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         ncwrite(save_file{iv},'rain',narr_z);
         ncwrite(save_file{iv},'rain_time',t);
         clear z_url narr_z;
         
      %-----------------------------------------------------------------%
      case 'lwrad'
            
         z1_url  = [base_url 'dlwrf.' num2str(year_on) '.nc']; 
         z2_url  = [base_url 'ulwrf.sfc.' num2str(year_on) '.nc']; 
         narr_z1 = ncread(z1_url,'dlwrf',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         narr_z2 = ncread(z2_url,'ulwrf',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         ncwrite(save_file{iv},'lwrad_down',narr_z1);
         ncwrite(save_file{iv},'lwrad_up',  narr_z2);
         ncwrite(save_file{iv},'lrf_time',t);
         clear z1_url z2_url narr_z1 narr_z2;
        
      %-----------------------------------------------------------------%
      case 'swrad'
            
         z1_url  = [base_url 'dswrf.' num2str(year_on) '.nc']; 
         z2_url  = [base_url 'uswrf.sfc.' num2str(year_on) '.nc']; 
         narr_z1 = ncread(z1_url,'dswrf',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         narr_z2 = ncread(z2_url,'uswrf',[narr_i0 narr_j0 i0_thismonth],[narr_ni narr_nj n_thismonth]); 
         ncwrite(save_file{iv},'swrad_down',narr_z1);
         ncwrite(save_file{iv},'swrad_up',  narr_z2);
         ncwrite(save_file{iv},'srf_time',t);
         clear z1_url z2_url narr_z1 narr_z2;
       
      %-----------------------------------------------------------------%   
      otherwise; error(['Unknown variable in var_to_get: ' var_to_get{iv}]);
          
    end  
    end
    
    % Onto next year and month
    month_on = month_on + 1;
    if(month_on==13) 
      year_on  = year_on+1; 
      month_on = 1;
    end

  end
       
end
%=========================================================================%
