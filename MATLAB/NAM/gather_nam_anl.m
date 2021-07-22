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
clc; addpath('../Utilities');   % clear mex; % (useful when debugging nc_gen_* codes)

% Options
date_start = [2017,01];     % Start year, month to gather data for
date_end   = [2019,01];     % End   year, month to gather data for (good to do +1 month for interpolation-in-time purposes)
max_wait   = 120;           % Max time to wait (secs) for web response (uses Java & parallelism, set to 0 to disable)
%data_dir   = '/gpfs/data/epscor/anelson5/OSOM_Data_Repo/NAM/';
data_dir   = 'F:/OSOM_Data_Repo/NAM/';  % My desktop
grid_file  = [data_dir 'nam_grid.nc'];  % Grid file
var_to_get = {'Uwind','Vwind','Pair','Tair','Qair','rain','lwrad_down','lwrad_up','swrad_down','swrad_up'};

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

% Information of possible variables to gather
% READ-IN NAME ; WRITE-OUT NAME ; # DIMS IN NAM FILE ; UNIT CONVERSION FACTOR ; TIME NAME
var_info = {'Uwind',        'u-component_of_wind_height_above_ground',  4,  '';         ...
            'Vwind',        'v-component_of_wind_height_above_ground',  4,  '';         ...
            'Pair',         'Pressure_reduced_to_MSL_msl',              3,  './100';    ...   
            'Tair',         'Temperature_height_above_ground',          4,  '-273.15';  ...
            'Qair',         'Relative_humidity_height_above_ground',    4,  '';         ...
            'rain',         'Total_precipitation_surface_*',            3,  '*';        ...     % *'s indicate a special case
            'lwrad_down',   'downward_long_wave_rad_flux_surface',      3,  '';         ...
            'lwrad_up',     'upward_long_wave_rad_flux_surface',        3,  '';         ...
            'swrad_down',   'downward_long_wave_rad_flux_surface',      3,  '';         ...
            'swrad_up',     'upward_long_wave_rad_flux_surface',        3,  '';         ...
            'Cfra',         'Total_cloud_cover_entire_atmosphere',      3,  '';         ...
            'sensible',     'Sensible_heat_net_flux_surface',           3,  '';         ...
            'latent',       'Latent_heat_net_flux_surface',             3,  '';         ...
           };
            
% Possible others to add: Tdew, albedo, surface drag coeff.
       
%=========================================================================%

% if grid file does not exist, ask user to do it
% if it does exist, load the grid variables
if(exist(grid_file,'file')~=2)
	error('Specified grid file does not exist.  Either fix the grid_file variable or run gather_nam_grid.m.');
else
	nam_mask = ncread(grid_file,'mask');
	nam_i0   = ncread(grid_file,'i0');
	nam_j0   = ncread(grid_file,'j0');
	[nam_ni,nam_nj] = size(nam_mask);
end

% ensure var_to_get values match possible ROMS surface forcing variables
n_var = numel(var_to_get);
i_var = zeros(n_var,1);
for iv=1:n_var
try  
    i_var(iv) = find( strcmp(var_info(:,1),var_to_get{iv}) == 1);
catch err
    error(['Variable in var_to_get ''' var_to_get{iv} ''' is not a valid ROMs surface forcing variable.  Check var_info for valid variables.']);
end
end
clear iv;
    
% Start year and month
year_on  = date_start(1);
month_on = date_start(2);
  
% Save directories
for iv=1:numel(var_to_get)
	var_dir = [data_dir var_to_get{iv}];
	if(exist(var_dir,'dir')~=7)
        mkdir(var_dir);
	end
end
clear iv var_dir;
  
%-------------------------------------------------------------------------%

% Loop through available months and years
while( (year_on + (month_on-0.5)/12) < (date_end(1) + (date_end(2)-0.25)/12) )

	% Number of days in the current month
	n_day = eomday(year_on,month_on);
    
    % Time vector
    t  = linspace(0,n_day-(1/8),n_day*8); % NAM data is in 3-hr files.  pre-2010, they are forecasts and don't include all variables.
    nt = numel(t);
    
    % Initialize save files for this year & month
    var_file = cell(n_var,1);
    for iv=1:n_var
    	var_file{iv} = [data_dir var_to_get{iv} '/NAM_' var_to_get{iv} '_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
        if(exist(var_file{iv},'file')~=2)
            nc_gen_nam_data(var_file{iv},nam_ni,nam_nj,{var_to_get{iv}});
        end
    end
    clear iv;
    
	% Loop through files for this month
	for it=1:nt
        
        % Time info
        day_on    = floor(t(it))+1;                 % valid values: 1:eomday(year_on,month_on)
        time_on   = floor( (t(it)-day_on+1)*4 )*6;  % valid values: 0, 6, 12, 18
        frcst_on  = (t(it)-day_on+1)*24 - time_on;	% valid values: 0, 3
    
        % Construct base url for this day
        this_url  = [];
        if(     datenum(year_on,month_on,day_on) >= datenum(date_min(1,:)) && datenum(year_on,month_on,day_on) <= datenum(date_max(1,:)) );     this_url = base_url{1};
        elseif( datenum(year_on,month_on,day_on) >= datenum(date_min(2,:)) && datenum(year_on,month_on,day_on) <= datenum(date_max(2,:)) );     this_url = base_url{2};
        else;   warning(['Date ' datestr(datenum(year_on,month_on,day_on)) ' does not have data, skipping...']);
        end

        % Continue if time frame is valid (url exists)
        if(~isempty(this_url))
        clc; disp(['Gathering NAM data for date ' datestr(datenum(year_on,month_on,day_on)) ' ' num2str(time_on) 'hrs...']);
        
            % Construct base URL string  
            this_url = [this_url num2str(year_on) sprintf('%0.2d',month_on) '/'];
            this_url = [this_url num2str(year_on) sprintf('%0.2d',month_on) sprintf('%0.2d',day_on) '/'];
            this_url = [this_url 'namanl_218_' num2str(year_on) sprintf('%0.2d',month_on) sprintf('%0.2d',day_on) '_' sprintf('%0.2d',time_on) '00_00'];
        
            % File URLs
            data_url1 = [this_url num2str(frcst_on)   '.grb'];
            data_url2 = [this_url num2str(frcst_on+3) '.grb'];	% occasionally forecast is needed needed for 'rain' variable
          
            % Info from these files
            info1 = ncinfo_web(max_wait,data_url1); info1_vars = {info1.Variables.Name};
            info2 = ncinfo_web(max_wait,data_url2); info2_vars = {info2.Variables.Name};
            
            % Loop through variables
            for iv=1:n_var
        
                % See if this variable at this time has already been gathered
                % If it hasn't, either time will be empty, or time at this index will be the fill value 
                check = false;
                time  = ncread(var_file{iv},'time');
                if(numel(time)<it)
                    check = true;   
                else
                    if(time(it)>1e36);  check = true;   end
                end
                
                % Continue if data hasn't been gathered
                if(check)
                switch var_to_get{iv}
                    
                    % Rain is a special case
                    case 'rain'
                        
                        % Rain variable always starts with this value
                        ii = find( contains(info2_vars,'Total_precipitation_surface_')==1 );
                        rain_name = info2_vars{ii};
                        
                        % Try to read variable
                        rain_value = ncread_web(max_wait, data_url2, rain_name, [nam_i0 nam_j0 1], [nam_ni nam_nj 1]);
                        
                        % Depending on rain name, different correction factor is needed
                        switch rain_name
                            case 'Total_precipitation_surface_3_Hour_Accumulation';             rain_value = rain_value./10800; % 3-hr accumulation into average rate per second
                            case 'Total_precipitation_surface_6_Hour_Accumulation';             rain_value = rain_value./21600; % 6-hr accumulation into average rate per second
                            case 'Total_precipitation_surface_Mixed_intervals_Accumulation';    rain_value = rain_value./10800;	% Slightly unsure about this one...
                            otherwise;                                                          error(['Unimplemented rain variable name:' rain_name]);
                        end
                        
                        % Save to file (& add time)
                        ncwrite(var_file{iv}, 'rain', rain_value, [1 1 it]);
                        ncwrite(var_file{iv}, 'time', t(it),      it      );
            
                        % Clean-up
                        clear ii rain_name rain_value 
                   
                    % All other variables are straightforward
                    otherwise
                        
                        % Make sure variable exists in online file before attempting to download it
                        if( any(strcmp(info1_vars,var_info{i_var(iv),2})==1) )
                            
                            % Variable dimensions in online file
                            var_dims1 = ones(var_info{i_var(iv),3},1);  var_dims1(1:2) = [nam_i0, nam_j0];
                            var_dims2 = ones(var_info{i_var(iv),3},1);  var_dims2(1:2) = [nam_ni, nam_nj];
                            
                            % Try to read variable
                            var_value = ncread_web(max_wait, data_url1, var_info{i_var(iv),2}, var_dims1, var_dims2);
                            
                            % Apply correction factor
                            eval(['var_value = var_value' var_info{i_var(iv),4} ';']);
                            
                            % Save to file (& add time)
                            ncwrite(var_file{iv}, var_to_get{iv}, var_value, [1 1 it]);
                            ncwrite(var_file{iv}, 'time',         t(it),     it      );
            
                            % Clean-up
                            clear var_dims1 var_dims2 var_value;

                        else
                            warning(['Variable ' var_info{i_var(iv),2} ' does not exist in file ' this_url ', skipping...']);
                        end

                    
                end
                end
                
                % Clean-up
                clear check time;
              
              
            end
            clear iv;
          
            % Clean-up
            clear this_url data_url data_url2
        
        end
        
        % Clean-up
        clear day_on time_on frcst_on;
        
	end
	clear it;
    
    % Clean-up
    clear n_day t nt var_file;

end
clear year_on month_on

%-------------------------------------------------------------------------%
if(0)   % OLD CODE

  % Loop through available months and years
  while( (year_on + (month_on-0.5)/12) < (date_end(1) + (date_end(2)-0.25)/12) )
    
    % Number of days for the current month
    n_day = eomday(year_on,month_on);
    
    % Time vector
    t = 0 : 0.125 : (n_day-1)+0.8751;
    nt = numel(t);
    
    % Initialize save files for this year & month
    save_file = cell(numel(var_to_get),1);
    for i=1:numel(var_to_get)
      save_file{i} = [data_dir var_to_get{i} '/NAM_' var_to_get{i} '_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
      if(exist(save_file{i},'file')~=2)
        nc_gen_nam_data(save_file{i},nam_ni,nam_nj,nt,{var_to_get{i}});
      end
    end
    clear i
    
%     % See if any data for this month has been read yet
%     test = ncread(save_file{end},'srf_time');
%     it0 = find(test>1e36,1,'first');
%     if(isempty(it0))
%       warning(['Data for date ' datestr(datenum(year_on,month_on,day_on)) ' has already been gathered, skipping...'])
%     else
    
      % Loop through files for this month
      for it=1:nt
        
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
          test = ncread(save_file{1},'wind_time',it,1);
          if(test>1e36)
            Uwind = ncread_web(max_wait,data_url,'u-component_of_wind_height_above_ground',[nam_i0 nam_j0 1 1],[nam_ni nam_nj 1 1]);
            Vwind = ncread_web(max_wait,data_url,'v-component_of_wind_height_above_ground',[nam_i0 nam_j0 1 1],[nam_ni nam_nj 1 1]);
            ncwrite(save_file{1},'Uwind',Uwind,[1 1 it]);   clear Uwind;
            ncwrite(save_file{1},'Vwind',Vwind,[1 1 it]);   clear Vwind;
            ncwrite(save_file{1},'wind_time',t(it),it);
          end
          
          % Pressure
          test = ncread(save_file{2},'pair_time',it,1);
          if(test>1e36)
            Pair = ncread_web(max_wait,data_url, 'Pressure_reduced_to_MSL_msl',           [nam_i0 nam_j0 1],  [nam_ni nam_nj 1]  ) ./ 100;
            ncwrite(save_file{2},'Pair',Pair,[1 1 it]);   clear Pair;
            ncwrite(save_file{2},'pair_time',t(it),it);
          end
            
          % Temeprature
          test = ncread(save_file{3},'tair_time',it,1);
          if(test>1e36)
            Tair = ncread_web(max_wait,data_url, 'Temperature_height_above_ground',       [nam_i0 nam_j0 1 1],[nam_ni nam_nj 1 1])  - 273.15;
            ncwrite(save_file{3},'Tair',Tair,[1 1 it]);   clear Tair;
            ncwrite(save_file{3},'tair_time',t(it),it);
          end
            
          % Dewpoint 
          test = ncread(save_file{4},'qair_time',it,1);
          if(test>1e36)  
            Qair = ncread_web(max_wait,data_url, 'Relative_humidity_height_above_ground', [nam_i0 nam_j0 1 1],[nam_ni nam_nj 1 1]);
            ncwrite(save_file{4},'Qair',Qair,[1 1 it]);   clear Qair;
            ncwrite(save_file{4},'qair_time',t(it),it);   
          end
          
          % Cloud fraction
          test = ncread(save_file{5},'cloud_time',it,1);
          if(test>1e36)  
            Cfra = ncread_web(max_wait,data_url, 'Total_cloud_cover_entire_atmosphere',   [nam_i0 nam_j0 1],  [nam_ni nam_nj 1]  );
            ncwrite(save_file{5},'Cfra',Cfra,[1 1 it]);   clear Cfra;
            ncwrite(save_file{5},'cloud_time',t(it),it);
          end
          
          % Rain (based on forecast)
          test = ncread(save_file{6},'rain_time',it,1);
          if(test>1e36) 
            try          rain = ncread_web(max_wait,data_url2,'Total_precipitation_surface_3_Hour_Accumulation', [nam_i0 nam_j0 1],  [nam_ni nam_nj 1]  ) ./ 10800;
            catch err
              try        rain = ncread_web(max_wait,data_url2,'Total_precipitation_surface_6_Hour_Accumulation', [nam_i0 nam_j0 1],  [nam_ni nam_nj 1]  ) ./ 21600; 
              catch err 
                  rain = ncread_web(max_wait,data_url2,'Total_precipitation_surface_Mixed_intervals_Accumulation', [nam_i0 nam_j0 1],  [nam_ni nam_nj 1]  ) ./ 10800; 
              end
            end
            ncwrite(save_file{6},'rain',rain,[1 1 it]);   clear rain;
            ncwrite(save_file{6},'rain_time',t(it),it);
          end
          
          % Longwave radiation
          test = ncread(save_file{7},'lrf_time',it,1);
          if(test>1e36) 
            lwrad_down = ncread_web(max_wait,data_url,'downward_long_wave_rad_flux_surface', [nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
            lwrad_up   = ncread_web(max_wait,data_url,'upward_long_wave_rad_flux_surface',   [nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
            ncwrite(save_file{7},'lwrad_down',lwrad_down,[1 1 it]);   clear lwrad_down;
            ncwrite(save_file{7},'lwrad_up',  lwrad_up,  [1 1 it]);   clear lwrad_up;
            ncwrite(save_file{7},'lrf_time',t(it),it);
          end
          
          % Shortwave radiation
          test = ncread(save_file{8},'srf_time',it,1);
          if(test>1e36) 
            swrad_down = ncread_web(max_wait,data_url,'downward_short_wave_rad_flux_surface',[nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
            swrad_up   = ncread_web(max_wait,data_url,'upward_short_wave_rad_flux_surface',  [nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
            ncwrite(save_file{8},'swrad_down',swrad_down,[1 1 it]);   clear swrad_down;
            ncwrite(save_file{8},'swrad_up',  swrad_up,  [1 1 it]);   clear swrad_up;
            ncwrite(save_file{8},'srf_time',t(it),it);
          end
          
        end
      end
      clear this_url data_url data_url2 day_on time_on;
        
    %end
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
