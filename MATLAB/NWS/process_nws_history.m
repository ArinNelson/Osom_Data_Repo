function process_nws_history(input_file,output_file)
%=========================================================================%
% process_nws_history(input_file,output_file)                             %
% Process data from raw NWS station data file, compute the daily          %
% statistics,  and write the results to the netcdf file output_file.      %
%                                                                         %
% process_nws_history(input_file,output_file,opt_name,opt_value,...)      %
% Specify any number of options:                                          %
%                                                                         %
% (More options will be added as needed.)                                 %
%                                                                         %
% Last updated by Arin Nelson on 09/20/2021                               %
%=========================================================================%

% % DEBUGGING
% clear mex;
% input_file  = './nws_TAN_raw.nc';
% output_file = './nws_TAN_daily.nc';

%-------------------------------------------------------------------------%
% INITIALIZATION

    % If output file not yet generated, generate it!
    if(exist(output_file,'file')~=2)
        %ncid   = netcdf.open(input_file,'NC_NOWRITE');
        %nws_id = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'station_id');
        nws_id = input_file(end-9:end-7);   % FOR NOW...
        gen_nws_file_stats(nws_id,output_file);
        %netcdf.close(ncid);
    end
    
    % Get index of last-written timestep
    ncid     = netcdf.open(output_file,'NC_NOWRITE');
    [~,i_on] = netcdf.inqDim(ncid,0);
    i_on     = i_on + 1;
    netcdf.close(ncid);
    
    % Load time variable
    time_datenum = ncread(input_file,'time');
    day_datenum  = floor(time_datenum);
    day_list     = day_datenum(1) : 1 : day_datenum(end);
    if(i_on <=numel(day_list))
        
        % Time info    
        i_start      = max([1 find(day_datenum>=day_list(i_on),1,'first')]);
        time_datenum = time_datenum(i_start:end);
        day_datenum  = day_datenum(i_start:end);
        day_list     = day_list(i_on:end);

        % Gather data from day_start onwards
        info      = ncinfo(input_file);
        var_name  = {info.Variables(5:end).Name};
        n_var     = numel(var_name);
        var_value = cell(n_var,1);
        for i=1:n_var
            var_value{i} = ncread(input_file,var_name{i},i_start,inf);
        end
        clear i info 
        
        % Loop through days
        for i = 1:numel(day_list)

%             % DEBUGGING
%             clc; disp(['On day ' num2str(i) '/' num2str(numel(day_list)) '...']);

            % Indices of data for this day
            i_day = find(day_datenum == day_list(i));
            if(~isempty(i_day))

                % Last data point from previous day & first data point from next day
                if(i_day(1)>1);         i_day = [find(day_datenum == day_list(i)-1,1,'last'); i_day(:)];   end
                if(i<numel(day_list)); i_day = [i_day(:); find(day_datenum == day_list(i)+1,1,'first')];	end

                % Loop through variables
                for iv=1:n_var
                
                    % This data
                    tmp_time  = time_datenum(i_day) - day_list(i);
                    tmp_value = var_value{iv}(i_day);
                    
                    % Sometimes 2+ samples are taken at the same time
                    [~,ii] = unique(tmp_time);
                    if(numel(ii)<numel(tmp_time))
                        jj = find(~ismember(1:numel(tmp_time),ii));
                        if(numel(jj)>0)
                            tmp_time(jj) = tmp_time(jj)+(1:numel(jj))'.*1e-9;
                        end
                    end
                    
                    % Continue if data is available
                    if(~all(isnan(tmp_value)))
                    
                        % Processing depends on variable
                        switch var_name{iv}

                            %-------------------------------------------------%
                            % Most variables are simple...
                            case {'tmpc','relh','mslp','cfra'}

                                % The count and percentiles only need data within the day
                                i0          = find(tmp_time>=0 & tmp_time<1);
                                tmp_num     = numel(i0);
                                tmp_prctile = prctile(tmp_value(i0),[0 10 25 50 75 90 100]);

%                                 % For now...
%                                 tmp_mean = nanmean(tmp_value(i0));
%                                 tmp_std  = nanstd(tmp_value(i0));
%                                 tmp_cvg  = tmp_time(i0(end)) - tmp_time(i0(1));
                                
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

                            %-------------------------------------------------%    
                            % Wind needs to have mean computed in complex plane
                            case {'sped'}

                                % Complex wind variable
                                % DEBUG CHECK: 
                                %   abs(wind_vec) - tmp_value < 1e-9;
                                %   (180/pi).*angle(wind_vec) - var_value{iv-1}(i_day) < 1e-9;
                                wind_vec = tmp_value.*exp(-sqrt(-1)*(pi/180).*var_value{iv-1}(i_day));

                                % The count and percentiles only need data within the day
                                i0          = find(tmp_time>=0 & tmp_time<1);
                                tmp_num     = numel(i0);
                                tmp_prctile = prctile(wind_vec(i0),[0 10 25 50 75 90 100]);

%                                 % For now...
%                                 tmp_mean = nanmean(wind_vec(i0));
%                                 tmp_std  = nanstd(wind_vec(i0));
%                                 tmp_cvg  = tmp_time(i0(end))-tmp_time(i0(1));
                                
                                % Mean and standard deviation are estimated from trapezoidal 
                                % method by taking  mean(flow) = trapz(flow)/secs_per_day
                                if(numel(tmp_time)>1)

                                    % Interp values at start and end of day
                                    ti = unique([0; tmp_time(tmp_time>0 & tmp_time<1); 1]);
                                    try
                                    xi = interp1(tmp_time,wind_vec,ti,'linear');
                                    catch err
                                       pause(1e-9); 
                                    end

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

                                % Save wind speed
                                if(exist('tmp_mean','var')==1)
                                if(~isempty(tmp_mean))    
                                    ncwrite(output_file,[var_name{iv} '_mean'],abs(tmp_mean),i_on);        
                                    ncwrite(output_file,[var_name{iv} '_std'], abs(tmp_std), i_on);        
                                    ncwrite(output_file,[var_name{iv} '_cvg'], tmp_cvg, i_on);              
                                end
                                end
                                ncwrite(output_file,[var_name{iv} '_num'],    tmp_num,i_on);                
                                ncwrite(output_file,[var_name{iv} '_prctile'],abs(tmp_prctile),[i_on 1]);	

                                % Save wind direction
                                if(exist('tmp_mean','var')==1)
                                if(~isempty(tmp_mean))    
                                    ncwrite(output_file,[var_name{iv-1} '_mean'],angle(tmp_mean).*(180/pi),i_on);       clear tmp_mean;
                                    ncwrite(output_file,[var_name{iv-1} '_std'], angle(tmp_std).*(180/pi));             clear tmp_std;
                                    ncwrite(output_file,[var_name{iv-1} '_cvg'], tmp_cvg, i_on);                        clear tmp_cvg;
                                end
                                end
                                ncwrite(output_file,[var_name{iv-1} '_num'],    tmp_num,i_on);                          clear tmp_num;
                                ncwrite(output_file,[var_name{iv-1} '_prctile'],angle(tmp_prctile).*(180/pi),[i_on 1]); clear tmp_prctile;

                            %-------------------------------------------------%    
                            % Rain is tricky since it's cumulative and resets at varying times
                            case {'p01m'}
                            if(~all(tmp_value(~isnan(tmp_value))==0))
                                tmp = tmp_value(~isnan(tmp_value));
                                ncwrite(output_file,[var_name{iv-1} '_mean'],   max(tmp),           i_on    );
                            	ncwrite(output_file,[var_name{iv-1} '_std'],    max(tmp)-min(tmp),	i_on    );
                                ncwrite(output_file,[var_name{iv-1} '_cvg'],    1,                  i_on    );
                                ncwrite(output_file,[var_name{iv-1} '_num'],    numel(tmp),         i_on    );                         
                                ncwrite(output_file,[var_name{iv-1} '_prctile'],prctile(tmp,[0 10 25 50 75 90 100]),	[i_on 1]);
                                clear tmp;
                            else
                                ncwrite(output_file,[var_name{iv-1} '_mean'],   0,         i_on    );
                            	ncwrite(output_file,[var_name{iv-1} '_std'],    0,         i_on    );
                                ncwrite(output_file,[var_name{iv-1} '_cvg'],    0,         i_on    );
                                ncwrite(output_file,[var_name{iv-1} '_num'],    0,         i_on    );                          clear tmp_num;
                                ncwrite(output_file,[var_name{iv-1} '_prctile'],zeros(1,7),[i_on 1]);
                            end
    
                        end
            
                    end
                    clear tmp_time tmp_value;
                    
                    % Save time
                    ncwrite(output_file,'time',day_list(i),i_on);
                    
                end
                clear iv;

            end
            clear i_day;

            % Onto next day
            i_on = i_on + 1;

        end
        clear i_day   

    end
    clear i_on;
        
end