function download_station_coops(save_file,coops_id,var_to_get)
%=========================================================================%
% download_station_coops(coops_id, save_file, var_to_get)                 %
% Download variables in cell array var_to_get from CO-OPS station and save%
% to NetCDF file save_file.                                               %
% by Arin Nelson                                                          %
% on 08/29/2021                                                           %
%                                                                         %
% var_to_get may include any of these, as long as they're in a cell array:%
%   water_level                                                           %
%   water_temperature                                                     %
%   salinity                                                              %
%   air_temperature                                                       %
%   air_pressure                                                          %
%   (others will be added as needed)                                      %
%                                                                         %
% download_station_coops(coops_id, save_file, option_name, ...            %
%                                             option_value, ...)          %
% Specify any number of options:                                          %
%                                                                         %
% 'Start Date'  (default: [1990 01 01])                                   %
%     Either a MATLAB datenum or a 1x3 vector specifying the              %
%     date as [year month day] to start looking for data from this gauge. %
%     Times are treated as being in the UTC timezone.                     %
%                                                                         %
% Last updated by Arin Nelson on 08/29/2021                               %
%=========================================================================%

%     % DEBUGGING
%     save_file = 'test.nc';
%     coops_id  = '8454000';
%     var_to_get = {'water_level','water_temperature','salinity','air_temperature','air_pressure'};

    % Max date to look for data for
    datenow = datetime(datestr(now),'Format','yyyy-MM-dd HH:mm z','TimeZone','local');

    % Validate required inputs
    if(~ischar(coops_id));      error('CO-OPS ID must be a string!');       end
    if(~ischar(save_file));     error('Save file name must be a string!');	end
    
    % Define option defaults
    start_date = [1990 01 01];
    
%     % Check if any optionsare specified
%     if(nargin>2)
%     for ia=1:2:nargin-3
%     switch varargin{ia}
%         %case lower('TS_IDs');           ts_id      = varargin{ia+1};
%         case lower('Start Date');       start_date = varargin{ia+1};
%     end
%     end
%     end
    
    % Save variables associated with input variables
    % (others will be added when needed)
    var_info = {'water_level',          'zeta', 'zeta_sigma'; ...
                'water_temperature',    'temp', ''; ... 
                'salinity',             'salt', 'spec_grav'; ...
                'air_temperature',      'tair', ''; ...
                'air_pressure',         'pair', ''; ...
               };
    
    % For reading csv from url as a table
    readtableweb = @(filename)readtable(filename);
    myoptions    = weboptions('ContentReader',readtableweb,'Timeout',600);
    
    %---------------------------------------------------------------------%

    % If save file doesn't yet exist, create it
    if(exist(save_file,'file')~=2)
        nc_gen_coops_stationfile(coops_id,save_file,var_to_get);
    end
    
    % Loop through variables
    n_var = numel(var_to_get);
    for iv=1:n_var
        
        % Index of variable in var_info
        ivar = find( strcmp(var_info(:,1),var_to_get{iv})==1 );
        
        % Get current time index
        old_time = ncread(save_file,[var_info{ivar,2} '_time']);
        if(~isempty(old_time))
            nn       = numel(old_time)+1;
            old_time = old_time(  end);
            yearI    = year( old_time); 
            monthI   = month(old_time);
            dayI     = day(  old_time);
        else
            nn       = 1;
            old_time = 0;
            yearI    = start_date(1);
            monthI   = start_date(2);
            dayI     = start_date(3);
        end
        
        % Final values
        timeF  = now;
        yearF  = year( timeF);
        monthF = month(timeF);
        dayF   = day(  timeF);  
        
        % Starts
        yearOn  = yearI;
        monthOn = monthI;
        dayOn   = dayI;
        
        % Loop through remaining years, months, and days
        while(yearOn <= yearF)
            
            % Last month/day
            if(yearOn < yearF)
                monthE = 12;
            else
                monthE = monthF;
                dayE   = dayF;
            end
            
            % Loop through months
            for im=monthOn:monthE
            
                % This timeframe's URL
                data_url = ['https://api.tidesandcurrents.noaa.gov/api/prod/datagetter' ...
                            '?begin_date=' num2str(yearOn) sprintf('%0.2d',im) sprintf('%0.2d',dayOn) ...
                            '&end_date='   num2str(yearOn) sprintf('%0.2d',im) sprintf('%0.2d',eomday(yearOn,im)) ...
                            '&station='    coops_id ...
                            '&datum=STND&time_zone=gmt&units=metric&format=csv&product=' var_info{ivar,1}];
            
                % Try to read the url
                try data_table = webread(data_url,myoptions);
                
                    % First column is date-time of measurement
                    this_time = datenum(data_table{:,1});
                    
                    % Only keep data beyond old_time
                    it = find(this_time>old_time,1,'first') : numel(this_time);
                    if(~isempty(it))
                        
                        % Truncate time
                        this_time = this_time(it);
                        
                        % Second column is data we want
                        this_data = [data_table{it,2}];
                        
                        % Some variables have a second data entry    
                        if(~isempty(var_info{iv,3}))
                            other_data = [data_table{:,3}];
                        end
                        
                        % Save variables
                        ncwrite(save_file,[var_info{ivar,2} '_time'],this_time,nn);
                        ncwrite(save_file,var_info{ivar,2},this_data,nn);
                        if(~isempty(var_info{iv,3}))
                            ncwrite(save_file,var_info{ivar,3},other_data,nn);
                        end

                        % Next timestep
                        nn = nn + numel(this_time);
                        
                        % Clean-up
                        clear this_data other_data;
                        
                    end
                    clear this_time it data_table;
                
                catch err    
                end   
                clear data_url;
                        
                % dayOn is now 1
                dayOn = 1;
                        
            end
            clear im monthE dayE
             
            % monthOn is now 1
            monthOn = 1;
                        
            % Next year
            yearOn  = yearOn+1;
            
        end
        clear ivar old_time yearOn monthOn dayOn timeF yearF monthF dayF;
        
    end
    clear iv;

end
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nc_gen_coops_stationfile(coops_id,file_name,var_name)

    % These URLs contain metadata for the station with the provided ID
    info_url  = ['https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' coops_id '.xml'];
    tides_url = ['https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' coops_id '/harcon.xml?units=metric'];
    
    % Station info
	info_raw        = xmlread_web(info_url);
	info_raw        = info_raw(1).Stations(1).Station(1);
	name            = info_raw.name.Text;
	lat             = str2double( info_raw.lat.Text );
	lon             = str2double( info_raw.lng.Text );
	affils          = info_raw.affiliations.Text;
	%is_tidal        = info_raw.tidal.Text;
	state           = info_raw.state.Text;
	tzone           = str2double( info_raw.timezonecorr.Text );
    
    % Tides info
    tides_raw = xmlread_web(tides_url);
    tides_raw = tides_raw(1).HarmonicConstituents(1).HarmonicConstituent;
	tide      = struct;
	if(numel(tides_raw)>1)
	for it=1:numel(tides_raw)
        tide(it).name      = tides_raw{it}.name.Text;
        tide(it).amplitude = str2double(tides_raw{it}.amplitude.Text);
        tide(it).phase_GMT = str2double(tides_raw{it}.phase_GMT.Text); 
    end
    end

    %---------------------------------------------------------------------%
    
    % Init the NetCDF
    ncid = netcdf.create(file_name,'NETCDF4');
    
    % NetCDF constants
    nc_global    = netcdf.getConstant('NC_GLOBAL');
    nc_unlimited = netcdf.getConstant('NC_UNLIMITED');
    
    % Some global attributes
    netcdf.putAtt(ncid,nc_global,'title',['NOAA CO-OPS Station Data, ID ' coops_id]);
    netcdf.putAtt(ncid,nc_global,'history',['Created on ' datestr(now)]);
    netcdf.putAtt(ncid,nc_global,'station_name',name);
    netcdf.putAtt(ncid,nc_global,'station_affiliations',affils);
    netcdf.putAtt(ncid,nc_global,'station_state',state);
    
    % Dimensions
    dimid = [];
    
    % Tide dimensions (KEEP EVEN IF EMPTY)
    dimid(end+1) = netcdf.defDim(ncid,'tide',numel(tide));
    dimid(end+1) = netcdf.defDim(ncid,'char',4);
    
    % Single-value variables
    varid = [];
    varid(end+1) = netcdf.defVar(ncid,'lat','double',[]);
        netcdf.putAtt(ncid,varid(end),'long_name','latitude');
        netcdf.putAtt(ncid,varid(end),'units','degrees_north');
    varid(end+1) = netcdf.defVar(ncid,'lon','double',[]);
        netcdf.putAtt(ncid,varid(end),'long_name','longitude');
        netcdf.putAtt(ncid,varid(end),'units','degrees_east');
    varid(end+1) = netcdf.defVar(ncid,'tzone','double',[]);
        netcdf.putAtt(ncid,varid(end),'long_name','timezone offset');
        netcdf.putAtt(ncid,varid(end),'units','hours from UTC');
        
    % Tide variables
    varid(end+1) = netcdf.defVar(ncid,'tide_name','char',[dimid(1) dimid(2)]);
        netcdf.putAtt(ncid,varid(end),'long_name','tidal constituent name');
        netcdf.putAtt(ncid,varid(end),'units','char');
    varid(end+1) = netcdf.defVar(ncid,'tide_amp','double',dimid(1));
        netcdf.putAtt(ncid,varid(end),'long_name','tidal constituent amplitude');
        netcdf.putAtt(ncid,varid(end),'units','meters above MSL');
    varid(end+1) = netcdf.defVar(ncid,'tide_phs','double',dimid(1));
        netcdf.putAtt(ncid,varid(end),'long_name','tidal constituent phase');
        netcdf.putAtt(ncid,varid(end),'units','degrees');
        netcdf.putAtt(ncid,varid(end),'relative_to','UTC');
        
    % Time variables
    for iv=1:numel(var_name)
    switch var_name{iv}
        case 'water_level';         dimid(end+1) = netcdf.defDim(ncid,'zeta_t',nc_unlimited);
        case 'water_temperature';   dimid(end+1) = netcdf.defDim(ncid,'temp_t',nc_unlimited);
        case 'salinity';            dimid(end+1) = netcdf.defDim(ncid,'salt_t',nc_unlimited);
        case 'air_temperature';     dimid(end+1) = netcdf.defDim(ncid,'tair_t',nc_unlimited);    
        case 'air_pressure';        dimid(end+1) = netcdf.defDim(ncid,'pair_t',nc_unlimited);  
    end
    end
    
    % Time-dependent variables
    %{'water_level','water_temperature','conductivity','salinity','currents',...
    %              'air_temperature','wind','air_pressure','humidity'};
    for iv=1:numel(var_name)
    switch var_name{iv} % MAKE SURE THESE ARE IN THE SAME ORDER AS THE TIME VARIABLES DEFINED ABOVE!
        case 'water_level'  % Includes 'Sigma' (uncertainty)
            varid(end+1) = netcdf.defVar(ncid,'zeta_time','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','water level time');
                netcdf.putAtt(ncid,varid(end),'units','matlab datenum');
                netcdf.putAtt(ncid,varid(end),'timezone','GMT');
            varid(end+1) = netcdf.defVar(ncid,'zeta','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','water level');
                netcdf.putAtt(ncid,varid(end),'units','meters above MSL');
                netcdf.putAtt(ncid,varid(end),'time','zeta_time');
            varid(end+1) = netcdf.defVar(ncid,'zeta_sigma','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','water level uncertainty');
                netcdf.putAtt(ncid,varid(end),'units','meters plus/minus from zeta');
                netcdf.putAtt(ncid,varid(end),'time','zeta_time');    
        case 'water_temperature'
            varid(end+1) = netcdf.defVar(ncid,'temp_time','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','water temperature time');
                netcdf.putAtt(ncid,varid(end),'units','matlab datenum');
                netcdf.putAtt(ncid,varid(end),'timezone','GMT');
            varid(end+1) = netcdf.defVar(ncid,'temp','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','water temperature');
                netcdf.putAtt(ncid,varid(end),'units','degrees_Celsius');
                netcdf.putAtt(ncid,varid(end),'time','temp_time'); 
        case 'salinity'     % Includes specific gravity
            varid(end+1) = netcdf.defVar(ncid,'salt_time','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','water salinity time');
                netcdf.putAtt(ncid,varid(end),'units','matlab datenum');
                netcdf.putAtt(ncid,varid(end),'timezone','GMT');
            varid(end+1) = netcdf.defVar(ncid,'salt','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','water salinity');
                netcdf.putAtt(ncid,varid(end),'units','psu');
                netcdf.putAtt(ncid,varid(end),'time','salt_time');
            varid(end+1) = netcdf.defVar(ncid,'spec_grav','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','water specific gravity');
                netcdf.putAtt(ncid,varid(end),'units','standard_units');
                netcdf.putAtt(ncid,varid(end),'time','salt_time');    
        case 'air_temperature'
            varid(end+1) = netcdf.defVar(ncid,'tair_time','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','air temperature time');
                netcdf.putAtt(ncid,varid(end),'units','matlab datenum');
                netcdf.putAtt(ncid,varid(end),'timezone','GMT');
            varid(end+1) = netcdf.defVar(ncid,'tair','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','air temperature');
                netcdf.putAtt(ncid,varid(end),'units','degrees_Celsius');
                netcdf.putAtt(ncid,varid(end),'time','tair_time');
        case 'air_pressure'
            varid(end+1) = netcdf.defVar(ncid,'pair_time','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','air pressure time');
                netcdf.putAtt(ncid,varid(end),'units','matlab datenum');
                netcdf.putAtt(ncid,varid(end),'timezone','GMT');
            varid(end+1) = netcdf.defVar(ncid,'pair','double',dimid(2+iv));
                netcdf.putAtt(ncid,varid(end),'long_name','air pressure');
                netcdf.putAtt(ncid,varid(end),'units','mb');
                netcdf.putAtt(ncid,varid(end),'time','tpair_time');
    end
    end
    
    % End definitions and close the netcdf file
    netcdf.endDef(ncid);
    netcdf.close(ncid);
    
    % Save available variables
    ncwrite(file_name,'lat',lat);
    ncwrite(file_name,'lon',lon);
    ncwrite(file_name,'tzone',tzone);
    if(numel(tides_raw)>1)
    	for it=1:numel(tide);  ncwrite(file_name,'tide_name',tide(it).name,[it 1]);   end
        ncwrite(file_name,'tide_amp',[tide.amplitude]);
        ncwrite(file_name,'tide_phs',[tide.phase_GMT]);
    end
    
end