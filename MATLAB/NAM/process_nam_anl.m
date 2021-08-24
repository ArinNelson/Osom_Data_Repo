%=========================================================================%
% process_nam_anl_v2.m
% Process North American Mesoscale (NAM) model Analysis (ANL) data
% by Arin Nelson
% on 07/18/2021
% 
% last edited by Arin Nelson on 08/14/2021
%=========================================================================%
clear; clc; addpath('../Utilities');   % clear mex;

% Time Options
date_start = [2018,01,01];	% Start year, month, day to gather data for
date_end   = [2018,12,31];	% End   year, month, day to gather data for

% Interpolation options
lon = -72.7 : 0.1 : -69.9;
lat =  40.5 : 0.1 :  42.2;

% Other options
save_file  = 'F:/OSOM_Data_Repo/ROMS/forcefiles/OSOM_frc_2018_NAM-ANL.nc';	% Save file name
save_title = 'NAM forcing, spatially variable';                         % Title attribute of frc file

% Shortwave radiation options
swrad_factor    = 1-0.23;             % Multiplicative factor (FROM DOPPIO)
swrad_dailyavg  = true;               % Save swrad as a daily average

% Variable options
var_to_get = {'wind','wind_time', '';       ...
              'Tair','tair_time', '';       ...
              'Pair','pair_time', '';       ...
              'Qair','qair_time', '';       ...
              'rain','rain_time', '';       ...
              'swrad','srf_time', 'net';    ...
              'lwrad','lrf_time', 'down';   ...
             };

% Source options
src_dir      = 'F:/OSOM_Data_Repo/NAM_ANL/';
src_gridfile = 'F:/OSOM_Data_Repo/NAM_ANL/namanl_grid.nc';
   
%-------------------------------------------------------------------------%

% Constants
src_dt = 3/24;
nx     = numel(lon);
ny     = numel(lat);

%=========================================================================%

% Interpolation meshgrid
[latm,lonm] = meshgrid(lat,lon);

% Data from grid file  
src_x    = ncread(src_gridfile,'x'   );
src_y    = ncread(src_gridfile,'y'   );
src_lon  = ncread(src_gridfile,'lon' );
src_lat  = ncread(src_gridfile,'lat' );
src_mask = ncread(src_gridfile,'mask');

% Source lon/lat meshgrid (if lon & lat are 1D)
if( any(size(src_lon)==1) )
	[tmpY, tmpX] = meshgrid(src_lat,src_lon);
	src_lon = tmpX;     clear tmpX;
	src_lat = tmpY;   	clear tmpY;
end
  
% Source x/y 2D grid (if x & y are 1D)
if( any(size(src_x)==1) )
    [tmpY, tmpX] = meshgrid(src_y,src_x);
	src_x = tmpX;      clear tmpX;
	src_y = tmpY;      clear tmpY;
end
  
% Estimate x and y at interpolation points
ntrplnt_x = scatteredInterpolant(src_lon(:),src_lat(:),src_x(:),'linear');
ntrplnt_y = scatteredInterpolant(src_lon(:),src_lat(:),src_y(:),'linear');
xi        = ntrplnt_x(lonm,latm);   clear ntrplnt_x;
yi        = ntrplnt_y(lonm,latm);   clear ntrplnt_y;

% Time array
src_t   = datenum(date_start) : src_dt : datenum(date_end)+1 - src_dt;
nperday = round(1/src_dt);

%-------------------------------------------------------------------------%

% Generate forcing file if it doesn't yet exist
n_var = size(var_to_get,1);
if(exist(save_file,'file')~=2)  
    tmp = var_to_get(:,1);
    for iv=1:n_var
    if( strcmp(var_to_get{iv,3},'down')==1 )
        tmp{iv} = [tmp{iv} '_down'];
    end
    end
	nc_gen_frc_roms(save_file,nx,ny,tmp,save_title);
	ncwrite(save_file,'lon',lon);
	ncwrite(save_file,'lat',lat);
    clear tmp;
end

% Dtermine timestep
t_test = ncread(save_file,var_to_get{end,2});
if(isempty(t_test))
    t_on  = datenum(date_start);
    nt_on = 1;
else
    t_on  = datenum(date_start) + src_dt.*numel(t_test);
    nt_on = numel(t_test)+1;
end

% Loop through years, months, days
imask = find(src_mask==1);
while(t_on <= datenum(date_end))
clc; disp(['On date ' datestr(t_on) '...']);

	% Year, month, day
	year_on  = num2str( year(t_on) );
	month_on = sprintf( '%0.2d', month(t_on) );
	day_on   = sprintf( '%0.2d', day(t_on) );
    
	% Loop through wanted variables
	for iv=1:n_var
            
        % Source data file
        src_file = [src_dir '/' var_to_get{iv} '/' year_on '/' month_on '/' var_to_get{iv} '_' year_on '_' month_on '_' day_on '.nc'];
            
        % Load data
        switch var_to_get{iv,1}
            
            case 'wind'
                src_data = ncread(src_file,'Uwind') + sqrt(-1).* ncread(src_file,'Vwind');
                
            case {'lwrad','swrad'}
                switch var_to_get{iv,3}
                    case 'net';  src_data = ncread(src_file,[var_to_get{iv,1} '_down']) - ncread(src_file,[var_to_get{iv,1} '_up']);
                    case 'down'; src_data = ncread(src_file,[var_to_get{iv,1} '_down']);
                    case 'up';   src_data = ncread(src_file,[var_to_get{iv,1} '_up']);
                end  
                
            otherwise
                src_data = ncread(src_file,var_to_get{iv,1});
    
        end
        src_time  = ncread(src_file,'time');
        
         % Depending on rain type, use a different correction factor
         if(strcmp(var_to_get{iv,1},'rain'))
             
            src_time = src_time(:,2:end);
            src_data = src_data(:,:,:,2:end);
            src_data(isnan(src_data)) = 0;
            src_type = squeeze(ncread(src_file,[var_to_get{iv,1} '_type'],[1 2 1],[inf inf inf]));
            for it=1:size(src_data,3)
            for ic=1:size(src_data,4)
                this_src = squeeze(src_type(it,ic,:));
                this_src(this_src>1e36) = [];
                this_src = char(this_src);
                    switch this_src'
                        case 'Total precipitation (3_Hour Accumulation) @ Ground or water surface'
                            src_data = src_data ./ (3*84600);   % 3-hr kg/m2 accumulation into kg/m2/s
                        case 'Total precipitation (Mixed_intervals Accumulation) @ Ground or water surface'
                            src_data = src_data ./ ( (src_time(it,ic)-src_time(it,ic-1))*84600);   % ?-hr kg/m2 accumulation into kg/m2/s
                        otherwise
                            pause(1);
                    end
                    clear this_src
            end
            end
            clear it ic src_type;
                
         end

        % Interpolate to output grid
        data = NaN(nx,ny,nperday);
        for id=1:nperday
            [it,ic] = find( src_time == round(src_dt*24)*(id-1) );
            if(numel(it)>1);    it=it(1);   ic=ic(1);   end
            if(~isempty(it))
                zz           = src_data(:,:,it,ic);  
                if( strcmp(var_to_get{iv,1},'wind') )
                    ntrplnt      = scatteredInterpolant(src_x(:),src_y(:),zz(:),'linear');
                else
                    ntrplnt      = scatteredInterpolant(src_x(imask),src_y(imask),zz(imask),'linear');
                end
                data(:,:,id) = ntrplnt(lonm,latm);
                clear zz ntrplnt;
            end
        end
        clear id it ic src_time;
        
        % Linearly interpolate data missing in time
        for ix=1:size(data,1)
        for iy=1:size(data,2)
            data(ix,iy,:) = nanfill(squeeze(data(ix,iy,:)));
        end
        end
        clear ix iy;
        
        % Some variables need a conversion factor
        if(strcmp(var_to_get{iv,1},'Tair'));    data = data - 273.15;   end     % Kelvin to Celsius
        if(strcmp(var_to_get{iv,1},'Pair'));    data = data ./ 100;     end     % Pa to mb
        if(strcmp(var_to_get{iv,1},'swrad'))
            if(~isempty(swrad_factor)); data = data.*swrad_factor;  end         % Apply DOPPIO correction factor
            if(swrad_dailyavg == true); data = nanmean(data,3);     end         % Compute daily average if specified  
        end
        
        % Quality check
        if(any(isnan(data(:))))
            pause(1);
        end
        
        % Save data
        switch var_to_get{iv,1}
            
            case 'wind'
                ncwrite(save_file,'Uwind',real(data),[1 1 nt_on]);
                ncwrite(save_file,'Vwind',imag(data),[1 1 nt_on]);
                
            case 'lwrad'
                switch var_to_get{iv,3}
                    case 'net';  src_data = ncread(src_file,'lwrad_down') - ncread(src_file,'lwrad_up');
                    case 'down'; src_data = ncread(src_file,'lwrad_down');
                    case 'up';   src_data = ncread(src_file,'lwrad_up'  );
                end  
                
            case 'swrad'
                if(swrad_dailyavg == true)
                switch var_to_get{iv,3}
                    case 'net';  ncwrite(save_file,'swrad',     data,[1 1 (nt_on-1)/(1/src_dt)+1]);
                    case 'down'; ncwrite(save_file,'swrad_down',data,[1 1 (nt_on-1)/(1/src_dt)+1]);
                    case 'up';   ncwrite(save_file,'swrad_up',  data,[1 1 (nt_on-1)/(1/src_dt)+1]);
                end
                else
                switch var_to_get{iv,3}
                    case 'net';  ncwrite(save_file,'swrad',     data,[1 1 nt_on]);
                    case 'down'; ncwrite(save_file,'swrad_down',data,[1 1 nt_on]);
                    case 'up';   ncwrite(save_file,'swrad_up',  data,[1 1 nt_on]);
                end
                end
                
            otherwise
                ncwrite(save_file,var_to_get{iv,1},data,[1 1 nt_on]);
    
        end

    end
    clear iv;
    
    % When done with variables, save time
    t_to_write = ((t_on) : src_dt : t_on+1-src_dt) - datenum(date_start);
    for iv=1:n_var
    if( swrad_dailyavg == true && strcmp(var_to_get{iv},'swrad')==1 )
        ncwrite(save_file,var_to_get{iv,2},t_to_write(1),(nt_on-1)/(1/src_dt)+1);
    else
        ncwrite(save_file,var_to_get{iv,2},t_to_write,nt_on);
    end
    end
    clear iv t_to_write;
    
    % Onto the next timestep
    t_on  = t_on  + 1;
    nt_on = nt_on + (1/src_dt);
        
    % Clean-up
    clear year_on month_on day_on
    
end
clear t_on 