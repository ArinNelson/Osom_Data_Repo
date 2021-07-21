%=========================================================================%
% gather_nam_anl.m
% Gather North American Mesoscale (NAM) model Analysis (ANL) data and sort
% by Arin Nelson
% on 07/15/2021
%
% Example file url for testing/debugging:
% dat_url = 'https://www.ncei.noaa.gov/thredds/dodsC/model-namanl-old/201001/20100101/namanl_218_20100101_0000_000.grb';
%
% Good reference:
% https://code.usgs.gov/coawstmodel/COAWST/-/blob/84738f369c84b51109c8f9a323c0ddbe9c6f7587/Tools/mfiles/mtools/ncei_2roms.m
% 
% last edited by Arin Nelson on 07/16/2021
%=========================================================================%
clc; clear mex; addpath('../Utilities');

% Switches
Switch    = zeros(9,1);
Switch(1) = 1;      % Initialize grid file and sampling info
Switch(2) = 1;      % Perform data gathering

% Options
date_start = [2017,01];     % Start year, month to gather data for
date_end   = [2021,01];     % End   year, month to gather data for
max_wait   = 120;           % Max time to wait (secs) for web response (set to 0 to disable)
data_dir   = '/gpfs/data/epscor/anelson5/OSOM_Data_Repo/NAM/';
var_incl   = {'winds','Pair','Tair','Qair','Cfra','rain','lwrad','swrad'};

%-------------------------------------------------------------------------%

% Constants
%date_min = [2004,03,02; 2020,05,18];   % Data before 2010 is not 3-hrly...
date_min = [2010,01,01; 2020,05,18];
date_max = [2020,05,15; 2020,12,31];
base_url = {'https://www.ncei.noaa.gov/thredds/dodsC/model-namanl-old/', ...
            'https://www.ncei.noaa.gov/thredds/dodsC/model-namanl/', ...
           };
grid_url = 'https://rda.ucar.edu/datasets/ds609.0/docs/latlon-g218.txt'; 
mask_url = 'https://www.ncei.noaa.gov/thredds/dodsC/model-namanl-old/201007/20100701/namanl_218_20100701_0000_000.grb';
lon_lim  = [-74.0 -69.0];
lat_lim  = [ 40.0  42.5];

% Information of possible variables to gather
% READ-IN NAME ; WRITE-OUT NAME ; # DIMS IN NAM FILE ; UNIT CONVERSION FACTOR ; TIME NAME
% var_info = {'u-component_of_wind_height_above_ground',          'Uair',         4,  '',         'wind_time';    ...
%             'v-component_of_wind_height_above_ground',          'Vair',         4,  '',         'wind_time';    ...
%             'Pressure_reduced_to_MSL_msl',                      'Pair',         3,  './100',    'pair_time';    ...     % PA to mbar
%             'Temperature_height_above_ground',                  'Tair',         4,  '-273.15',  'tair_time';    ...     % K to C
%             'Relative_humidity_height_above_ground',            'Qair',         4,  '',         'qair_time';    ...
%             'Total_cloud_cover_entire_atmosphere',              'Cfra',         3,  '',         'cloud_time';   ...
%             'Total_precipitation_surface_3_Hour_Accumulation',  'rain',         3,  './10800',  'rain_time';    ...     % NOTE: THIS IS FROM 3HR FORECAST!
%             'downward_short_wave_rad_flux_surface',             'swrad_down',   3,  '',         'srf_time';     ...
%             'upward_short_wave_rad_flux_surface',               'swrad_up',     3,  '',         'srf_time';     ...
%             'downward_long_wave_rad_flux_surface',              'lwrad_down',   3,  '',         'lrf_time';     ...
%             'upward_long_wave_rad_flux_surface',                'lwrad_up',     3,  '',         'lrf_time';     ...
%            };
 
%=========================================================================%
if(Switch(1))

  % Delete grid file if it already exists  
  if(exist(grid_file,'file')==2)
    delete(grid_file);
  end
    
  % Read in NAM 218 lon/lat grid (from url)
  tmp1 = textscan( webread(grid_url,weboptions('ContentType','text')) , '%s', 'delimiter', '\n' );
  tmp1 = tmp1{1};
  tmp1(1:2) = [];
  tmp2 = zeros(4,numel(tmp1));
  for i=1:numel(tmp1)
    tmp3 = textscan(tmp1{i},'%f %f %f %f');
    tmp2(:,i) = [tmp3{:}];
  end
  clear tmp1 tmp3 i;
  
  % Write as 2D grid
  nx      = max(tmp2(1,:));
  ny      = max(tmp2(2,:));
  nam_lon = zeros(nx,ny);
  nam_lat = zeros(nx,ny);
  for i=1:size(tmp2,2)
    nam_lat(tmp2(1,i),tmp2(2,i)) =  tmp2(3,i);
    nam_lon(tmp2(1,i),tmp2(2,i)) = -tmp2(4,i); % Units are in deg.W, so -deg.E
  end
  clear tmp2 i nx ny;  
    
  % Find grid over specified region
  [ii,jj] = find(nam_lon>=lon_lim(1) & nam_lon<=lon_lim(end) & nam_lat>lat_lim(1) & nam_lat<lat_lim(end));
  ii      = min(ii)-1 : max(ii)+1;
  jj      = min(jj)-1 : max(jj)+1;
  
  % These variables are what will be used by ncread, etc.
  nam_i0 = ii(1);  ni = numel(ii);
  nam_j0 = jj(1);  nj = numel(jj);
    
  % Truncate grid
  nam_lon = nam_lon(ii,jj);
  nam_lat = nam_lat(ii,jj);
  
  % Read in additional grid variables
  nam_x    = ncread_web(max_wait,mask_url,'x',nam_j0,ni);
  nam_y    = ncread_web(max_wait,mask_url,'y',nam_j0,nj);
  nam_mask = ncread_web(max_wait,mask_url,'Land_cover_land1_sea0_surface',[nam_i0 nam_j0 1],[ni nj 1]);
    
  % Generate grid file
  nc_gen_nam_grid(grid_file,ni,nj);
  
  % Save these variables to it
  ncwrite(grid_file,'x',   nam_x   );
  ncwrite(grid_file,'y',   nam_y   );
  ncwrite(grid_file,'lon', nam_lon );
  ncwrite(grid_file,'lat', nam_lat );
  ncwrite(grid_file,'mask',nam_mask);
  ncwrite(grid_file,'i0',  nam_i0  );
  ncwrite(grid_file,'j0',  nam_j0  );
  
end
%=========================================================================%
if(Switch(2))
    
  % Load necessary variables
  nam_mask = ncread(grid_file,'mask');
  nam_i0   = ncread(grid_file,'i0');
  nam_j0   = ncread(grid_file,'j0');
  [nam_ni,nam_nj] = size(nam_mask);
    
  % Start year and month
  year_on  = date_start(1);
  month_on = date_start(2);
  
  % Save directories
  for i=1:numel(var_incl)
    this_dir = [data_dir var_incl{i}];
    if(exist(this_dir,'dir')~=7)
      mkdir(this_dir);
    end
  end
  
  % Loop through available months and years
  while( (year_on + (month_on-0.5)/12) < (date_end(1) + (date_end(2)-0.25)/12) )
    
    % Number of days for the current month
    n_day = eomday(year_on,month_on);
    
    % Time vector
    t = 0 : 0.125 : (n_day-1)+0.8751;
    nt = numel(t);
    
    % Initialize save files for this year & month
    save_file = cell(numel(var_incl),1);
    for i=1:numel(var_incl)
      save_file{i} = [data_dir var_incl{i} '/NAM_' var_incl{i} '_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
      if(exist(save_file{i},'file')~=2)
        nc_gen_nam_data(save_file{i},nam_ni,nam_nj,nt,{var_incl{i}});
      end
    end
    clear i
    
    % See if any data for this month has been read yet
    test = ncread(save_file{end},'srf_time');
    it0 = find(test>1e36,1,'first');
    if(isempty(it0))
      warning(['Data for date ' datestr(datenum(year_on,month_on,day_on)) ' has already been gathered, skipping...'])
    else
    
      % Loop through files for this month
      for it=it0:nt
        
        % Time info
        day_on    = floor(t(it))+1;
        time_on   = floor( (t(it)-day_on+1)*4 )*6;
        frcst_on  = (t(it)-day_on+1)*24 - time_on;
        
        % Construct base url for this day
        this_url  = [];
        if( datenum(year_on,month_on,day_on) >= datenum(date_min(1,1),date_min(1,2),date_min(1,3)) & ...
            datenum(year_on,month_on,day_on) <= datenum(date_max(1,1),date_max(1,2),date_max(1,3)) )
          this_url = base_url{1};
        elseif( datenum(year_on,month_on,day_on) >= datenum(date_min(2,1),date_min(2,2),date_min(2,3)) & ...
                datenum(year_on,month_on,day_on) <= datenum(date_max(2,1),date_max(2,2),date_max(2,3)) )
          this_url = base_url{2};
        else
          warning(['Date ' datestr(datenum(year_on,month_on,day_on)) ' does not have data, skipping...']);
        end
      
        % Continue if time frame is valid
        if(~isempty(this_url))
        clc; disp(['Gathering NAM data for date ' datestr(datenum(year_on,month_on,day_on)) '...']);
        
          % Construct base URL string  
          this_url = [this_url num2str(year_on) sprintf('%0.2d',month_on) '/'];
          this_url = [this_url num2str(year_on) sprintf('%0.2d',month_on) sprintf('%0.2d',day_on) '/'];
          this_url = [this_url 'namanl_218_' num2str(year_on) sprintf('%0.2d',month_on) sprintf('%0.2d',day_on) '_' sprintf('%0.2d',time_on) '00_00'];
        
          % File URLs
          data_url  = [this_url num2str(frcst_on)   '.grb'];
          data_url2 = [this_url num2str(frcst_on+3) '.grb']; % for precip
          
          % Ensure files exist before continuing
          %test_avail = false;
          %try info=ncinfo(data_url); test_avail = true; catch err; end
          %if(test_avail==true)
              
          % Download winds
          Uair = ncread_web(max_wait,data_url,'u-component_of_wind_height_above_ground',[nam_i0 nam_j0 1 1],[nam_ni nam_nj 1 1]);
          Vair = ncread_web(max_wait,data_url,'v-component_of_wind_height_above_ground',[nam_i0 nam_j0 1 1],[nam_ni nam_nj 1 1]);
        
          % Meteorology
          Pair = ncread_web(max_wait,data_url, 'Pressure_reduced_to_MSL_msl',           [nam_i0 nam_j0 1],  [nam_ni nam_nj 1]  ) ./ 100;
          Tair = ncread_web(max_wait,data_url, 'Temperature_height_above_ground',       [nam_i0 nam_j0 1 1],[nam_ni nam_nj 1 1])  - 273.15;
          Qair = ncread_web(max_wait,data_url, 'Relative_humidity_height_above_ground', [nam_i0 nam_j0 1 1],[nam_ni nam_nj 1 1]);
          Cfra = ncread_web(max_wait,data_url, 'Total_cloud_cover_entire_atmosphere',   [nam_i0 nam_j0 1],  [nam_ni nam_nj 1]  );
          
          % Rain (based on forecast)
          try
            rain = ncread_web(max_wait,data_url2,'Total_precipitation_surface_3_Hour_Accumulation', [nam_i0 nam_j0 1],  [nam_ni nam_nj 1]  ) ./ 10800;
          catch err
            rain = ncread_web(max_wait,data_url2,'Total_precipitation_surface_Mixed_intervals_Accumulation', [nam_i0 nam_j0 1],  [nam_ni nam_nj 1]  ) ./ 10800;  
          end
          
          % Radiations
          lwrad_down = ncread_web(max_wait,data_url,'downward_long_wave_rad_flux_surface', [nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
          lwrad_up   = ncread_web(max_wait,data_url,'upward_long_wave_rad_flux_surface',   [nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
          swrad_down = ncread_web(max_wait,data_url,'downward_short_wave_rad_flux_surface',[nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
          swrad_up   = ncread_web(max_wait,data_url,'upward_short_wave_rad_flux_surface',  [nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
        
          % Save data
          ncwrite(save_file{1},'wind_time',t(it),it);
          ncwrite(save_file{1},'Uair',Uair,[1 1 it]);   clear Uair;
          ncwrite(save_file{1},'Vair',Vair,[1 1 it]);   clear Vair;
          ncwrite(save_file{2},'pair_time',t(it),it);
          ncwrite(save_file{2},'Pair',Pair,[1 1 it]);   clear Pair;
          ncwrite(save_file{3},'tair_time',t(it),it);
          ncwrite(save_file{3},'Tair',Tair,[1 1 it]);   clear Tair;
          ncwrite(save_file{4},'qair_time',t(it),it);   
          ncwrite(save_file{4},'Qair',Qair,[1 1 it]);   clear Qair;
          ncwrite(save_file{5},'cloud_time',t(it),it);
          ncwrite(save_file{5},'Cfra',Cfra,[1 1 it]);   clear Cfra;
          ncwrite(save_file{6},'rain_time',t(it),it);
          ncwrite(save_file{6},'rain',rain,[1 1 it]);   clear rain;
          ncwrite(save_file{7},'lrf_time',t(it),it);
          ncwrite(save_file{7},'lwrad_down',lwrad_down,[1 1 it]);   clear lwrad_down;
          ncwrite(save_file{7},'lwrad_up',  lwrad_up,  [1 1 it]);   clear lwrad_up;
          ncwrite(save_file{8},'srf_time',t(it),it);
          ncwrite(save_file{8},'swrad_down',swrad_down,[1 1 it]);   clear swrad_down;
          ncwrite(save_file{8},'swrad_up',  swrad_up,  [1 1 it]);   clear swrad_up;
          
          %end
          %clear test_avail
          
        end
      end
      clear this_url data_url data_url2 day_on time_on;
        
    end
    clear it nt t save_file test it0;
      
    % Onto next year and month
    month_on = month_on + 1;
    if(month_on==13) 
      year_on  = year_on+1; 
      month_on = 1;
    end
      
  end
  clear year_on month_on;
    
end
%=========================================================================%
