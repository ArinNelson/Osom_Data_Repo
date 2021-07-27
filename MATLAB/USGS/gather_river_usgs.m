function status = gather_river_usgs(usgs_id,year_range,save_dir)
% status = gather_river_usgs(usgs_id, date_range, save_dir)
% Gather river gauge data from specified USGS river gauge for the provided
% year range and save it to the specified directory.
% 
% Inputs
%   usgs_id:    USGS 8-digit ID number given as a string
%   year_range: Array of integer year numbers (e.g., 2010:2020)
%   save_dir:   Directory where downloaded data will be written
% 
% by Arin Nelson
% on 07/20/2021
% Last edited 07/27/2021
%=========================================================================%

    % FOR DEBUGGING
    %usgs_id    = '01118500';
    %year_range = 2018:2020;
    %save_dir   = 'D:/OSOM_Data_Repo/USGS/RiverGauges/';

    % URL for station info 
    info_url = ['https://waterservices.usgs.gov/nwis/site/?format=rdb&sites=' usgs_id];
    
    % URL for station inventory
    inv_url  = ['https://waterdata.usgs.gov/nwis/inventory/?site_no=' usgs_id '&format=rdb'];

    % Get data (also ensures it's a valid USGS station)
    try
        [info_name, info_value] = webread_usgs(info_url);
        [inv_name,  inv_value]  = webread_usgs(inv_url);
    catch err
        error(['Either given station ID ' usgs_id ' does not exist or the web service is down.']);
    end
    
    % Init save directory
    save_dir = [save_dir '/' usgs_id '/'];
    if(exist(save_dir,'dir')~=7); mkdir(save_dir); end
    
    % Save info file if it doesn't yet exist
    info_file = [save_dir usgs_id '_info.mat'];
    if(exist(info_file,'file')~=2)
        
        % Download station info
        station_name  = info_value{ strcmp(info_name,'station_nm'   )==1 };    % Station name
        longitude     = info_value{ strcmp(info_name,'dec_long_va'  )==1 };    % Decimal longitude
        latitude      = info_value{ strcmp(info_name,'dec_lat_va'   )==1 };    % Decimal latitude
        altitude      = info_value{ strcmp(info_name,'alt_va'       )==1 };    % Altitude in feet above NAV88
        drainage_area = inv_value{  strcmp(inv_name, 'drain_area_va')==1 };    % Drainage area in square miles
        
        
        % Save data
        save(info_file,'station_name','longitude','latitude','altitude','drainage_area');

    end
    
    % Loop through requested years and get list of available variables and the
	% the total number of samples available for each year
    n_year = numel(year_range);
    for iy=1:n_year
        
        % Construct data url for this station and year
        data_url = ['https://nwis.waterservices.usgs.gov/nwis/iv/?format=rdb' ...
                    '&sites=' usgs_id ...
                    '&startDT=' num2str(year_range(iy)) '-01-01' ...
                    '&endDT=' num2str(year_range(iy)) '-12-31' ...
                    '&siteStatus=all'];
        
        % Read in complete data record for this year at this station       
        [var_name, var_value] = webread_usgs(data_url);   
   
        % Gather time
        i_time   = find( strcmp( var_name, 'datetime'  )==1 );  
        all_time = datenum(var_value(:,i_time))';
        
        % Remove time from full data array
        var_name(   i_time) = [];
        var_value(:,i_time) = [];
        
        % Save the variable names if they're not empty
        if(~all(cellfun(@isempty,var_value)))
       
            % Init save dir
            this_dir = [save_dir num2str(year_range(iy)) '\'];
            if(exist(this_dir,'dir')~=7);  mkdir(this_dir);    end

            % Gather data from each variable for this year
            for iv=1:numel(var_name)
            if(exist([this_dir var_name{iv} '.mat'],'file')~=2) 
                
                % List of values
                this_list = var_value(:,iv);
                
                % Set empty entries to NaN
                i_empty = find( cellfun(@isempty,this_list)==1 );
                if(~isempty(i_empty))
                    for ie=1:numel(i_empty);   this_list{i_empty(ie)} = NaN;    end
                end
                
                % Variable values
                this_value = [this_list{:}];
                
                % Keep time & data when its available
                ii    = find(~isnan(this_value));
                time  = all_time(ii);
                value = this_value(ii);
                
                % Save to file
                save([this_dir var_name{iv} '.mat'],'time','value');
                
                % Clean-up
                clear this_list i_empty ie this_value ii time value;
                
            end
            end
            clear iv;
            
        end
        clear data_url var_name var_value all_time;
        
    end
    clear iy;

end