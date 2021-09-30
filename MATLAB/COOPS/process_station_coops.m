function process_station_coops(data_dir,coops_id)

%     % FOR DEBUGGING
%     data_dir = 'F:/OSOM_Data_Repo/Stations/COOPS/';
%     coops_id = '8454000';

    % Variable info (load name, save name)
    var_info = {'pair',         'pair_time';    ...
                'tair',         'tair_time';    ...
                'salt',         'salt_time';    ...
                'spec_grav',    'salt_time';	...
                'zeta',         'zeta_time';    ...
                'temp',         'temp_time';    ...
               };
                
    % Gather station info and available variables
    file_name = [data_dir '/COOPS_' coops_id '.nc'];
    
    % Parse out wanted info
    lat   = ncread(file_name,'lat');
    lon   = ncread(file_name,'lon');
    tzone = ncread(file_name,'tzone');
    
    % Eval string
    save_str = ['save(''' data_dir '/COOPS_' coops_id '_stats.mat'',''lon'',''lat'',''tzone'''];
    
    % Loop through variables and collect their timeseries
    n_var = size(var_info,1);
    for iv=1:n_var
        eval([var_info{iv,1} '=ncread(file_name,''' var_info{iv,1} ''');']);
        eval([var_info{iv,2} '=ncread(file_name,''' var_info{iv,2} ''');']);
    end
    clear iv;
    
    % Get min and max day from each time variable
    tLim = [inf 0];
    for iv=1:n_var
        eval(['tmp_time='  var_info{iv,2} ';']);
        tLim(1) = min([tLim(1) tmp_time(1)]);
        tLim(2) = max([tLim(2) tmp_time(end)]);
    end
    clear iv tmp_time;
    
    % Generate hourly time series
    stats_time = (floor(tLim(1)) : (1/24) : ceil(tLim(2))-(1/24)) + (1/48);
    nt         = numel(stats_time);
    save_str   = [save_str  ',''stats_time'''];
    tt         = floor(stats_time*24);
    
    % Compute hourly statistics of each variable
    for iv=1:n_var
        
        % Gather this variable
        eval(['tmp_time=floor('  var_info{iv,2} '.*24);']);
        eval(['tmp_value=' var_info{iv,1} ';']);
        
        % Init stats variable
        tmp_stats = NaN(nt,6);
        
        % Loop through hours
        for it=1:nt
        clc; disp(['On Variable ' num2str(iv) '/' num2str(n_var) ', hour ' num2str(it) '/' num2str(nt) '...']);
        
            % Indices for this hour
            ii = find(tmp_time>=tt(it)-(1/48),1,'first') : find(tmp_time<tt(it)+(1/48),1,'last');
        
            % If there's data available, compute the statistics
            if(~isempty(ii))
                vv = tmp_value(ii);   
                vv = vv(~isnan(vv));
                if(~isempty(vv))
                    tmp_stats(it,1) = numel(vv);
                    tmp_stats(it,2) = min(vv);
                    tmp_stats(it,3) = max(vv);
                    tmp_stats(it,4) = median(vv);
                    tmp_stats(it,5) = mean(vv);
                    tmp_stats(it,6) = std(vv);
                end
            end
        end
        save_str = [save_str ',''' var_info{iv,1} '_stats'''];
        eval([var_info{iv,1} '_stats = tmp_stats;']);
        clear tmp_stats tmp_time tmp_value it ii vv;
    end
    clear iv;
    
    % Save data to master file
    save_str = [save_str ');'];
    eval(save_str);

end