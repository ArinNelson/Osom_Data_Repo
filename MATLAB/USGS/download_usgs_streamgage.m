function download_usgs_streamgage(usgs_id,save_file,varargin)
%=========================================================================%
% download_usgs_streamgage(usgs_id, save_file)                            %
% Download data from USGS river gauge and save to NetCDF file save_file.  %
%                                                                         %
% download_usgs_streamgage(usgs_id, save_file, opt_name, opt_value, ...)  %
% Specify any number of options:                                          %
%                                                                         %
% (Options will be added as needed)                                       %
%                                                                         %
% Last updated by Arin Nelson on 09/16/2021                               %
%=========================================================================%

%     % DEBUG
%     usgs_id   = '01108410';
%     save_file = ['./usgs_' usgs_id '_raw.nc'];

    % Max date to look for data for
    datenow = datetime(datestr(now),'Format','yyyy-MM-dd HH:mm z','TimeZone','local');

    % Validate required inputs
    if(~ischar(usgs_id));	error('USGS ID must be a string!');         end
    if(~ischar(save_file));	error('Save file name must be a string!');	end

%     % Check if any options are specified
%     if(nargin>2)
%     for ia=1:2:nargin-3
%     switch varargin{ia}
% 
%     end
%     end
%     end
    
    %---------------------------------------------------------------------%
    
    % If save file doesn't yet exist, create it
    if(exist(save_file,'file')~=2)
        if(nargin>2)
            gen_usgs_file_raw(usgs_id,save_file,varargin{1});
        else
            gen_usgs_file_raw(usgs_id,save_file);
        end
    end
    
    % Get latest available date
    time_old = ncread(save_file,'time');
    if(isempty(time_old))
        info     = ncinfo(save_file);
        tmp      = info.Attributes( strcmp({info.Attributes.Name},'established date')==1 ).Value;
        tmp      = [tmp(1:4) '-' tmp(5:6) '-' tmp(7:8) ' 00:00:00'];
        last_day = datetime(tmp,'TimeZone','UTC');
        n_on     = 1;
    else
        last_day = datetime(datestr(time_old(find(time_old>0,1,'last'))+1e-9),'TimeZone','UTC');
        n_on     = numel(time_old);
    end
    
    % Loop through times up to current day
    for iy = year(last_day) : 1 : year(datenow)
     
%         % DEBUGGING
%         clc; disp(['On year ' num2str(year(last_day)) '...']);
        
        % Construct data URL
        data_url = ['https://nwis.waterservices.usgs.gov/nwis/iv/?format=rdb' ...
                    '&sites=' usgs_id ...
                    '&startDT=' datestr(last_day-1,'yyyy-mm-dd') ...
                    '&endDT='   num2str(year(last_day)+1) '-01-01' ...
                    '&siteStatus=all'];
        
        % Download data from URL as a table
        var_table = webread_usgs(data_url);
        var_names = var_table.Properties.VariableNames;
        if(~all(cellfun(@isempty,var_table.Variables)))
        
            % Get time in terms of matlab datetime
            this_time = NaT(size(var_table,1),1,'Format','yyyy-MM-dd HH:mm z','TimeZone','UTC');
            for ie=1:size(var_table,1)
                this_time(ie) = datetime([var_table.datetime{ie} ' ' var_table.tz_cd{ie}],'Format','yyyy-MM-dd HH:mm z','TimeZone','UTC');
            end
            clear ie;
            
            % Remove duplicate dates
            [this_time,ii] = unique(this_time);
            var_table = var_table(ii,:);
            clear ii;
            
        	% Ensure to only keep data for this_time >= date_on
            i_get     = find(this_time > last_day,1,'first') : numel(this_time);
            var_table = var_table(i_get,:);
            this_time = this_time(i_get);
            
            % Gather data
            n_entries = size(var_table,1);
            if(n_entries>0)
                
                % Save variables
                for iv=5:2:numel(var_names)
                    this_var = var_names{iv}(end-4:end);
                    switch this_var
                        case '00010';   this_save = 'temperature';
                        case '00060';   this_save = 'discharge';
                        case '00065';   this_save = 'gage_height';
                        case '00095';   this_save = 'conductance';
                        case '00300';   this_save = 'dissolved_oxygen';
                        case '00400';   this_save = 'potential_hydrogen';
                    end
                    ncwrite(save_file,this_save,cellfun(@str2double,var_table{:,iv}),n_on);
                    this_char = char(ones(n_entries,1));
                    for i=1:n_entries
                    if(~isempty(var_table{i,iv+1}{1}))    
                        this_char(i) = var_table{i,iv+1}{1}(1);
                    end
                    end
                    ncwrite(save_file,[this_save '_qcf'],this_char,n_on);
                end
            
                % Save time
                ncwrite(save_file,'time',datenum(this_time),n_on);
            
                % On to the next year
                n_on     = n_on + numel(this_time);
                last_day = datetime(datestr(datenum(this_time(end))),'TimeZone','UTC');
                
            else
                last_day = datetime(datestr(datenum([year(last_day)+1,1,1])),'TimeZone','UTC');
            end
            
            % Clean-up
            clear i_get iv this_time this_var i this_char

        else
            last_day = datetime(datestr(datenum([year(last_day)+1,1,1])),'TimeZone','UTC');
        end
        clear var_table data_url;
    
    end
    clear iy;

end