function download_nws_history(nws_id,save_file,date_begin,varargin)
%=========================================================================%
% download_nws_history(nws_id, save_file)                                 %
% Download data from NWS weather station and save to netcdf file          %
% save_file from start_date (as MATLAB datenum) up to the most recent     %
% value.  Data from archive at https://mesonet.agron.iastate.edu          %
%                                                                         %
% download_nws_history(nws_id, save_file, opt_name, opt_value, ...)       %
% Specify any number of options:                                          %
%                                                                         %
% (Options will be added as needed)                                       %
%                                                                         %
% Last updated by Arin Nelson on 09/19/2021                               %
%=========================================================================%

%     % DEBUGGING
%     nws_id     = 'TAN';
%     date_begin = '1990-01-01';
%     save_file  = ['./nws_' nws_id '_raw.nc'];
    
    % Max date to look for data for
    date_now = datetime(datestr(now),'Format','yyyy-MM-dd HH:mm z','TimeZone','local');

    % Validate required inputs
    if(~ischar(nws_id));	error('USGS ID must be a string!');         end
    if(~ischar(save_file));	error('Save file name must be a string!');	end
    
    %---------------------------------------------------------------------%
    
    % If save file doesn't yet exist, create it
    if(exist(save_file,'file')~=2)
        gen_nws_file_raw(nws_id,save_file);
    end
    
    % Get latest saved date
    time_old = ncread(save_file,'time');
    if(isempty(time_old) || all(isnan(time_old)))
        last_day = datetime(date_begin,'TimeZone','UTC');
        n_on     = 1;
    else
        last_day = datetime(datestr(time_old(find(time_old>0,1,'last'))-1e-9),'TimeZone','UTC');
        n_on     = numel(time_old);
    end
    
    % Init start and end date
    start_date = last_day;
    end_date   = date_now;
    
    % Loop
    date_start = start_date;
    date_end   = datetime(datestr(min([datenum(year(start_date),12,31),datenum(end_date)])),'TimeZone','UTC');
    var_to_get = {'tmpc','relh','drct','sped','mslp','p01m','skyc1','skyc2','skyc3'};
    while(year(date_start)<=year(date_end))
    
%         % DEBUGGING
%         clc; disp(['On year ' num2str(year(date_start)) '...']);
        
        % Construct data URL
        nws_url = ['https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=' nws_id];
        for iv=1:numel(var_to_get);             nws_url = [nws_url '&data=' var_to_get{iv}];	end
        if(any(strcmp(var_to_get,'skyc1')==1)); nws_url = [nws_url '&data=metar'];              end
        nws_url = [nws_url '&year1=' num2str(year(date_start)) '&month1=' num2str(month(date_start)) '&day1=' num2str(day(date_start))	...
                           '&year2=' num2str(year(date_end  )) '&month2=' num2str(month(date_end  )) '&day2=' num2str(day(date_end  ))	...
                    '&latlon=no&elev=no&missing=empty&trace=0.0001&direct=no&report_type=1&report_type=2'];  
        
        % Load data within specified time span
        readtableweb = @(filename)readtable(filename);
        nws_table    = webread(nws_url,weboptions('ContentType','text','Timeout',600,'ContentReader',readtableweb));
        
        % If now empty, there's no data available
        if(isempty(nws_table))
            %warning(['No data exists for station ' nws_id ' for the specified time frame, skipping...']);
        else
        
            % Get available variables
            nws_vars = nws_table.Properties.VariableNames;

            % Get time variable
            nws_time = datenum(nws_table.valid);
            
            % Keep data after last_day
            ii        = find(nws_time >= datenum(date_start));
            nws_time  = nws_time(ii);
            nws_table = nws_table(ii,:);
            nt        = numel(nws_time);
            
            % Save variables to master file
            for iv=1:numel(nws_vars)
            switch nws_vars{iv}
            
            	% Most variables are straightforward to transfer
                case {'tmpc','mslp','relh','p01m','drct','sped'}    
                    eval(['tmp = nws_table.' nws_vars{iv} ';']);
                    ncwrite(save_file,nws_vars{iv},tmp,n_on);  
                    clear tmp;
                    
                % Cloud fraction requires some more work...
                case 'skyc1'
                    if(~iscell(nws_table.skyc1));  tmp1 = cell(numel(nws_table.skyc1),1); else;   tmp1=nws_table.skyc1;  end
                    if(~iscell(nws_table.skyc2));  tmp2 = cell(numel(nws_table.skyc2),1); else;   tmp2=nws_table.skyc2;  end
                    if(~iscell(nws_table.skyc3));  tmp3 = cell(numel(nws_table.skyc3),1); else;   tmp3=nws_table.skyc3;  end
                    cstr = [tmp1, tmp2, tmp3];
                    Mstr = nws_table.metar;
                    cfra = NaN(nt,1);
                    for it=1:nt
                    if(~all(isempty([cstr{it,:}])))

                       % Available cloud values
                       tmp = NaN(3,1);
                       for i=1:3
                       if(~isempty(cstr{it,i}))    
                       switch cstr{it,i}
                           case 'CLR';  tmp(i) = 0;
                           case 'FEW';  tmp(i) = 1.5/8;
                           case 'SCT';  tmp(i) = 3.5/8;
                           case 'BKN';  tmp(i) = 5.5/8;
                           case 'OVC';  tmp(i) = 1;
                       end
                       end
                       end
                       
                       % Take max value
                       cfra(it) = nanmax(tmp);
                       clear tmp;
                       
                    end
                    end
                    
                    % Save
                    ncwrite(save_file,'cfra',cfra,n_on);
                    clear tmp1 tmp2 tmp3 cstr Mstr cfra it i;
                    
                % Others?
                otherwise
                    % sorry nothing
            end
            end
            clear iv;
            
%             % Debug look?
%             if(1)
%                 subplot(3,2,1); plot(nws_time,nws_table.tmpc,'.b'); title('tmpc');
%                 subplot(3,2,2); plot(nws_time,nws_table.relh,'.b'); title('relh');
%                 subplot(3,2,3); plot(nws_time,nws_table.drct,'.b'); title('drct');
%                 subplot(3,2,4); plot(nws_time,nws_table.sped,'.b'); title('sped');
%                 subplot(3,2,5); plot(nws_time,nws_table.mslp,'.b'); title('mslp');
%                 subplot(3,2,6); plot(nws_time,nws_table.p01m,'.b'); title('p01m');
%             end
            
            % Save time variable
            ncwrite(save_file,'time',nws_time,n_on);
            n_on = n_on + nt;
            
        end
            
        % Onto next year
        date_start = datetime(datestr(datenum(year(date_start)+1,1,1)),'TimeZone','UTC');
        date_end   = datetime(datestr(min([datenum(year(date_start),12,31),datenum(end_date)])),'TimeZone','UTC');
        
    end
    
end