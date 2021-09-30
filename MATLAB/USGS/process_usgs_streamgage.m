function process_usgs_streamgage(input_file,output_file)
%=========================================================================%
% process_usgs_streamgage_daily(input_file,output_file)                   %
% Process data from raw streamgage data file, compute the daily           %
% statistics, and write the results to the netcdf file output_file.       %
%                                                                         %
% process_usgs_streamgage(input_file,output_file,opt_name,opt_value,...)  %
% Specify any number of options:                                          %
%                                                                         %
% (More options will be added as needed.)                                 %
%                                                                         %
% Last updated by Arin Nelson on 09/16/2021                               %
%=========================================================================%

% % DEBUGGING
% clear mex;
% input_file  = './usgs_01108410_raw.nc';
% output_file = './usgs_01108410_daily.nc';

%-------------------------------------------------------------------------%
% INITIALIZATION

    % If output file not yet generated, generate it!
    if(exist(output_file,'file')~=2)
        nc_glbl = netcdf.getConstant('NC_GLOBAL');
        ncid    = netcdf.open(input_file,'NC_NOWRITE');
        usgs_id = netcdf.getAtt(ncid,nc_glbl,'site identification number');
        gen_usgs_file_stats(usgs_id,output_file);
        netcdf.close(ncid);
    end
    
    % Get index of last-written timestep
    ncid     = netcdf.open(output_file,'NC_NOWRITE');
    [~,i_on] = netcdf.inqDim(ncid,0);
    i_on     = i_on + 1;
    netcdf.close(ncid);
    
    % Load time variable
    time_raw  = ncread(input_file,'time');
    time_day  = floor(time_raw);
    time_days = min(time_day) : 1 : max(time_day);
    if(i_on <=numel(time_days))
        
        % Time info    
        i_start   = max([1 find(time_day>=time_days(i_on),1,'first')-1]);
        time_raw  = time_raw(i_start:end);
        time_day  = time_day(i_start:end);
        time_days = time_days(i_on:end);

        % Gather data from day_start onwards
        info      = ncinfo(input_file);
        var_name  = {info.Variables(6:2:end).Name};
        n_var     = numel(var_name);
        var_value = cell(n_var,1);
        var_qcf   = cell(n_var,1);
        for i=1:n_var

            % Read in data
            var_value{i} = ncread(input_file,var_name{i},         i_start,inf);
            var_qcf{i}   = ncread(input_file,[var_name{i} '_qcf'],i_start,inf);

            % Set data with poor qcf to NaN
            var_value{i}( ~ismember(var_qcf{i},'A') ) = NaN;

        end
        clear i info 

        % Loop through days
        for i = 1:numel(time_days)

    %         % DEBUGGING
    %         clc; disp(['On day ' num2str(i) '/' num2str(numel(time_days)) '...']);

            % Indices of data for this day
            i_day = find(time_day == time_days(i));
            if(~isempty(i_day))

                % Last data point from previous day & first data point from next day
                if(i_day(1)>1);         i_day = [find(time_day == time_days(i)-1); i_day(:)];   end
                if(i<numel(time_days)); i_day = [i_day(:); find(time_day == time_days(i)+1)];	end

                % Loop through variables
                for iv=1:n_var

                    % This data
                    tmp_time  = time_raw(i_day) - time_days(i);
                    tmp_value = var_value{iv}(i_day);

                    % Remove NaN data
                    tmp_time( isnan(tmp_value)) = [];
                    tmp_value(isnan(tmp_value)) = [];

                    % The count and percentiles only need data within the day
                    i0          = find(tmp_time>=0 & tmp_time<1);
                    tmp_num     = numel(i0);
                    tmp_prctile = prctile(tmp_value(i0),[0 10 25 50 75 90 100]);

                    % Mean and standard deviation are estimated from trapezoidal 
                    % method by taking  mean(flow) = trapz(flow)/secs_per_day
                    if(numel(tmp_time)>1)

                        % Interp values at start and end of day
                        ti = unique([0; tmp_time(tmp_time>0 & tmp_time<1); 1]);
                        xi = interp1(tmp_time,tmp_value,ti,'linear');

                        % Compute
                        if(sum(~isnan(xi))>1)
                            if(any(isnan(xi))); xi(isnan(xi)) = interp1(ti(~isnan(xi)),xi(~isnan(xi)),ti(isnan(xi)),'linear');  end
                            ti = ti(~isnan(xi));
                            xi = xi(~isnan(xi));
                            [tmp_mean, tmp_std] = mytrapz(ti,xi);
                            tmp_cvg = ti(end)-ti(1);
                        elseif(sum(~isnan(xi))==1)
                            tmp_mean = xi(~isnan(xi));
                            tmp_std  = inf;
                            tmp_cvg  = 1e-9;
                        end

                    else
                        tmp_mean = tmp_value;
                        tmp_std  = inf;
                        tmp_cvg  = 1e-9;
                    end

                    % Save
                    if(exist('tmp_mean','var')==1)
                    if(~isempty(tmp_mean))    
                        ncwrite(output_file,[var_name{iv} '_mean'],tmp_mean,i_on);          clear tmp_mean;
                        ncwrite(output_file,[var_name{iv} '_std'], tmp_std, i_on);          clear tmp_std;
                        ncwrite(output_file,[var_name{iv} '_cvg'], tmp_cvg, i_on);          clear tmp_cvg;
                    end
                    end
                    ncwrite(output_file,[var_name{iv} '_num'],    tmp_num,i_on);            clear tmp_num;
                    ncwrite(output_file,[var_name{iv} '_prctile'],tmp_prctile,[i_on 1]);	clear tmp_prctile;
                    ncwrite(output_file,'time',time_days(i),i_on);

                end
                clear iv;

            end

            % Onto next day
            i_on = i_on + 1;

        end
        clear i_day
    
    end
    clear i_on

end