%=========================================================================%
% gather_nam_anl.m
% Gather North American Mesoscale (NAM) model Analysis (ANL) data and sort
% by Arin Nelson
% on 07/15/2021
%
% Example file url for testing/debugging:
% dat_url = 'https://www.ncei.noaa.gov/thredds/dodsC/model-namanl-old/201001/20100101/namanl_218_20100101_0000_000.grb';
% 
% NOTE: Files from different time periods have different formats, variable names, etc...
% oldest file:      https://www.ncei.noaa.gov/thredds/dodsC/model-namanl-old/200403/20040302/namanl_218_20040302_1800_000.grb;
%   - variables related to lwrad, swrad, u, v are all lower-case
% last 'old' file:      https://www.ncei.noaa.gov/thredds/dodsC/model-namanl-old/202005/20200515/namanl_218_20200515_0000_000.grb2
%   - at some point, switch from grb to grb2. This is why reading in the list of available files first is recommended.
%   - at some point radiations become capitalized, but u and v are still lower-case
% first 'modern' file:  https://www.ncei.noaa.gov/thredds/dodsC/model-namanl/202005/20200518/nam_218_20200518_0000_000.grb2
%   - similar to last 'old' file, includes forecasts for +1hr, +2hrs
% last 'modern' file: same
% 
% NOTE: This NAM product has max 3-hrly resolution, the 218 product has hrly resolution but doesn't begin until 2019/06/18.
% 
% Good reference from ROMS variables:
% https://code.usgs.gov/coawstmodel/COAWST/-/blob/84738f369c84b51109c8f9a323c0ddbe9c6f7587/Tools/mfiles/mtools/ncei_2roms.m
% 
% last edited by Arin Nelson on 08/14/2021
%=========================================================================%
clear; clc; addpath('../Utilities');   % clear mex; % (useful when debugging nc_gen_* codes)

% Options
date_start = [2021,07,12];     % Start year, month, day to gather data for
date_end   = [2021,08,01];     % End   year, month, day to gather data for
max_wait   = 600;              % Max time to wait (secs) for web response (uses Java & parallelism, set to 0 to disable)
%data_dir   = '/gpfs/data/epscor/anelson5/OSOM_Data_Repo/NAM_ANL/';
data_dir   = 'F:/OSOM_Data_Repo/NAM_ANL/';  % My desktop
grid_file  = [data_dir 'namanl_grid.nc'];  % Grid file
var_to_get = {'lwrad','swrad','wind','Pair','Tair','Qair','rain','Cfra'};

%-------------------------------------------------------------------------%

% Constants
date_min = [2004,03,02; 2020,05,18];
date_max = [2020,05,15; 2021,09,01];
base_url = 'https://www.ncei.noaa.gov/thredds/';
data_url = {'model-namanl-old/', ...
            'model-namanl/', ...
           };

% Information of possible variables to gather
% READ-IN NAME ; WRITE-OUT NAME ; # DIMS IN NAM FILE ; UNIT CONVERSION FACTOR ; TIME NAME
var_info = {'wind'          'wind_time',    '-component_of_wind_height_above_ground',           4,  '';         ...
            'Pair',         'pair_time',    'Pressure_reduced_to_MSL_msl',                      3,  './100';    ...   
            'Tair',         'tair_time',    'Temperature_height_above_ground',                  4,  '-273.15';  ...
            'Qair',         'qair_time',    'Relative_humidity_height_above_ground',            4,  '';         ...
            'rain',         'rain_time',    'Total_precipitation_surface_',                     3,  '*';        ...     % *'s indicate a special case
            'lwrad',        'lrf_time',     '_Long-Wave_Radp_Flux_surface',                     3,  '';         ...
            'swrad',        'srf_time',     '_Short-Wave_Radiation_Flux_surface',               3,  '';         ...
            'Cfra',         'cloud_time',   'Total_cloud_cover_entire_atmosphere_single_layer',	3,  '';         ...
            'sensible',     'sen_time',     'Sensible_heat_net_flux_surface',                   3,  '';         ...
            'latent',       'lat_time',     'Latent_heat_net_flux_surface',                     3,  '';         ...
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

% Save directories
for iv=1:numel(var_to_get)
	var_dir = [data_dir var_to_get{iv}];
	if(exist(var_dir,'dir')~=7)
        mkdir(var_dir);
	end
end
clear iv var_dir;

%-------------------------------------------------------------------------%

% Start year, month, day
day_start = datenum(date_start);
day_end   = datenum(date_end  );

% Initialzie time variable
t_on  = day_start;

% Loop through available months and years
while( t_on <= day_end )
clc; disp(['Gathering data for date ' datestr(t_on) ]);  

    % Time variables
    year_on  = num2str(year(t_on));
    month_on = sprintf('%0.2d',month(t_on));
    day_on   = sprintf('%0.2d',day(t_on));
    
    % Loop through variables
    for iv=1:n_var
    
        % Check for existance of save directories, and make them if need be
        save_dir = [data_dir '/' var_to_get{iv} '/' year_on];
        if(exist(save_dir,'dir')~=7);   mkdir(save_dir);    end
        save_dir = [save_dir '/' month_on];
        if(exist(save_dir,'dir')~=7);   mkdir(save_dir);    end
        
        % Check for existance of save file
        save_file = [save_dir '/' var_to_get{iv} '_' year_on '_' month_on '_' day_on '.nc'];
        if(exist(save_file,'file')~=2)
            
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . %
            % If not yet done so, read URL for list of files available for this day
            if(exist('file_list','var')~=1)
            disp(' Gathering data file info...');    
                
                % Construct catalog URL
                catalog_url = [base_url 'catalog/'];
                if(     datenum(t_on) >= datenum(date_min(1,:)) && datenum(t_on) <= datenum(date_max(1,:)) );   catalog_url = [catalog_url data_url{1}];
                elseif( datenum(t_on) >= datenum(date_min(2,:)) && datenum(t_on) <= datenum(date_max(2,:)) );   catalog_url = [catalog_url data_url{2}];
                else;   warning(['Date ' datestr(t_on) ' does not have data, skipping...']);
                end
                catalog_url = [catalog_url year_on month_on '/' year_on month_on day_on '/catalog.xml'];
                
                % Read files available at this url
                try
                  file_list = ls_tds(catalog_url);
                  file_list = file_list(end:-1:1,:);  % Files are given in reverse order.  Fix that.
                  n_files   = size(file_list,1);
                catch err
                end
                if(exist('file_list','var')==1)
                    
                    % Gather available analysis times and forecast times
                    file_info = cell(n_files,4);
                    for i=1:n_files
                        tmp = strsplit(file_list{i,2},{'_','.'});                                   % Split file name into pieces at _ and .
                        file_info{i,1} = str2double(tmp{end-2}(1:2));                                                % Analysis time
                        file_info{i,2} = str2double(tmp{end-1}(2:3));                                              % Forecast time
                        file_info{i,3} = tmp{end};                                                  % File extension
                    end
                    clear i;

                    % File info #4 is the file's info structure
                    ii = [];
                    for i=1:n_files
                        this_url       = [base_url 'dodsC/' file_list{i,1}];
                        try
                        file_info{i,4} = ncinfo_web(max_wait,this_url);
                        catch err; ii = [ii, i];
                        end     
                    end
                    if(~isempty(ii));   file_list(ii,:) = []; file_info(ii,:) = []; n_files = size(file_info,1); end

                    % Gather analysis times and forecast times based on each
                    avail_main = unique([file_info{:,1}]);
                    n_main     = numel(avail_main);
                    avail_cast = cell(n_main,1);
                    n_cast     = zeros(n_main,1);
                    for i=1:n_main
                        ii = find( [file_info{:,1}]==avail_main(i) );
                        avail_cast{i} = [file_info{ii,2}];
                        n_cast(i)     = numel(ii);
                    end
                    clear i ii;

                    % The time variable will be 2D
                    nt_main = n_main;
                    nt_cast = max(n_cast);
                
                end

            end
            if(exist('file_list','var')==1)
                
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . %
                % Depending on variable name, gather necessary data
                disp([' On ' var_to_get{iv} '...']);

                % Construct name(s) of variable(s) to read in
                data_var_name = var_info{i_var(iv),3};

                % Init data values
                switch var_to_get{iv}
                    case {'wind','lwrad','swrad'}
                        data = NaN(nam_ni,nam_nj,nt_main,nt_cast,2);
                    case 'rain'
                        data = NaN(nam_ni,nam_nj,nt_main,nt_cast);
                        rain_type = cell(nt_main,nt_cast); 
                    otherwise
                        data = NaN(nam_ni,nam_nj,nt_main,nt_cast);
                end
                time_value = NaN(nt_main,nt_cast);

                % Loop through available files
                for i=1:n_files

                    % See if data is available in this file
                    switch var_to_get{iv}
                        case 'lwrad'
                            try0 = {file_info{i,4}.Variables.Name};
                            try1 = '_long_wave_rad_flux_surface';
                            try2 = '_Long-Wave_Radp_Flux_surface'; 
                            try3 = '_Long-Wave_Radiation_Flux_surface'; 
                            ii   = find( contains(lower(try0),lower(try1))==1 | contains(lower(try0),lower(try2))==1 | contains(lower(try0),lower(try3))==1 );
                        case 'swrad'
                            try0 = {file_info{i,4}.Variables.Name};
                            try1 = '_short_wave_rad_flux_surface';
                            try2 = '_Short-Wave_Radp_Flux_surface'; 
                            try3 = '_Short-Wave_Radiation_Flux_surface';
                            ii   = find( contains(lower(try0),lower(try1))==1 | contains(lower(try0),lower(try2))==1 | contains(lower(try0),lower(try3))==1 );
                        otherwise
                            ii = find( contains( lower({file_info{i,4}.Variables.Name}), lower(data_var_name) )==1 );
                    end

                    % If data available, get it!
                    if(~isempty(ii))

                        % Construct data file URL
                        this_url = [base_url 'dodsC/' file_list{i,1}];

                        % Data time index
                        imain = find( file_info{i,1} == avail_main );
                        icast = find( file_info{i,2} == avail_cast{imain} );
                        time_value(imain,icast) = avail_main(imain) + avail_cast{imain}(icast);

                        % Read in data
                        switch var_to_get{iv}
                            case 'wind'
                                jj = find( contains( lower({file_info{i,4}.Variables(ii).Name}), 'u-' )==1 );
                                data(:,:,imain,icast,1) = ncread(this_url,file_info{i,4}.Variables(ii(jj)).Name,[nam_i0 nam_j0 1 1],[nam_ni nam_nj 1 1]);
                                jj = find( contains( lower({file_info{i,4}.Variables(ii).Name}), 'v-' )==1 );
                                data(:,:,imain,icast,2) = ncread(this_url,file_info{i,4}.Variables(ii(jj)).Name,[nam_i0 nam_j0 1 1],[nam_ni nam_nj 1 1]);
                            case {'lwrad','swrad'}
                                jj = find( contains( lower({file_info{i,4}.Variables(ii).Name}), 'downward' )==1 );
                                data(:,:,imain,icast,1) = ncread(this_url,file_info{i,4}.Variables(ii(jj)).Name,[nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
                                jj = find( contains( lower({file_info{i,4}.Variables(ii).Name}), 'upward' )==1 );
                                data(:,:,imain,icast,2) = ncread(this_url,file_info{i,4}.Variables(ii(jj)).Name,[nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
                            case 'rain'
                                data(:,:,imain,icast)  = ncread(this_url,file_info{i,4}.Variables(ii).Name,[nam_i0 nam_j0 1],[nam_ni nam_nj 1]);
                                rain_type{imain,icast} = file_info{i,4}.Variables(ii).Attributes(1).Value;
                            case 'Tair'
                                sstr = '';  for j=1:numel(ii);  sstr(end+1) = file_info{i,4}.Variables(ii(j)).Name(1); end
                                jj = find(sstr=='T');
                                ddim = ones(var_info{i_var(iv),4},1);   ddim(1) = nam_i0;   ddim(2) = nam_j0;
                                ndim = ones(var_info{i_var(iv),4},1);   ndim(1) = nam_ni;   ndim(2) = nam_nj;
                                data(:,:,imain,icast) = ncread(this_url,file_info{i,4}.Variables(ii(j)).Name,ddim,ndim);
                            case 'Qair'
                                sstr = '';  for j=1:numel(ii);  sstr(end+1) = file_info{i,4}.Variables(ii(j)).Name(1); end
                                jj = find(sstr=='R');
                                ddim = ones(var_info{i_var(iv),4},1);   ddim(1) = nam_i0;   ddim(2) = nam_j0;
                                ndim = ones(var_info{i_var(iv),4},1);   ndim(1) = nam_ni;   ndim(2) = nam_nj;
                                data(:,:,imain,icast) = ncread(this_url,file_info{i,4}.Variables(ii(j)).Name,ddim,ndim);
                            otherwise
                                ddim = ones(var_info{i_var(iv),4},1);   ddim(1) = nam_i0;   ddim(2) = nam_j0;
                                ndim = ones(var_info{i_var(iv),4},1);   ndim(1) = nam_ni;   ndim(2) = nam_nj;
                                data(:,:,imain,icast) = ncread(this_url,file_info{i,4}.Variables(ii).Name,ddim,ndim);
                        end

                    end
                    clear ii;


                end
                clear i;

                % After data gathered, create the file and save the data
                switch var_to_get{iv}
                    case 'wind'
                        nc_gen_nam_data(save_file,nam_ni,nam_nj,nt_main,nt_cast,var_to_get(iv));
                        ncwrite(save_file,'Uwind',data(:,:,:,:,1));
                        ncwrite(save_file,'Vwind',data(:,:,:,:,2));
                    case 'lwrad'   
                        nc_gen_nam_data(save_file,nam_ni,nam_nj,nt_main,nt_cast,var_to_get(iv));
                        ncwrite(save_file,'lwrad_down',data(:,:,:,:,1));
                        ncwrite(save_file,'lwrad_up',  data(:,:,:,:,2));
                    case 'swrad'   
                        nc_gen_nam_data(save_file,nam_ni,nam_nj,nt_main,nt_cast,var_to_get(iv));
                        ncwrite(save_file,'swrad_down',data(:,:,:,:,1));
                        ncwrite(save_file,'swrad_up',  data(:,:,:,:,2));
                    case 'rain'
                        nc_gen_nam_data(save_file,nam_ni,nam_nj,nt_main,nt_cast,var_to_get(iv));
                        ncwrite(save_file,'rain',data);
                        for ii=1:nt_main
                        for jj=1:nt_cast
                        if(~isempty(rain_type{ii,jj}))
                            tmp = zeros(1,1,length(rain_type{ii,jj}));
                            tmp(1,1,:) = char(rain_type{ii,jj});
                            ncwrite(save_file,'rain_type',tmp,[ii jj 1]);
                        end
                        end
                        end
                        clear rain_type ii jj tmp;
                    otherwise      
                        nc_gen_nam_data(save_file,nam_ni,nam_nj,nt_main,nt_cast,var_to_get(iv));
                        ncwrite(save_file,var_to_get{iv},data);
                end

                % Save time variable
                ncwrite(save_file,'time',time_value);

                % Clean-up
                clear tdata time_value;

            end
        end
        
        % Clean-up
        clear save_dir save_file data_var_name;
        
    end
    clear iv;
    
    % Onto the next day
    t_on = t_on + 1;
    
    % Clean-up
    clear year_on month_on day_on file_list;
    
end
clear day_start day_end day_on;