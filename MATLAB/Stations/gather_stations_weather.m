% gather_kpvd.m
% Gather obs from weather stations in/around OSOM domain
% by Arin Nelson
% on 08/01/2021
% 
% Data provided by https://mesonet.agron.iastate.edu/request/download.phtml
% 
% last edited on 08/01/2021
%=========================================================================%
clc; addpath('../Utilities'); clear;

% Options
%wthr_to_get = {'BID','UUU','OQU','SVZ','PVD','WST'};
%date_min    = [1990,01,01];
%date_max    = [2020,12,31];
year_to_get = 2018:2020;
save_dir    = 'F:/OSOM_Data_Repo/Stations/Weather/';
var_to_get  = {'tmpc','relh','drct','sped','mslp','p01m','skyc1','skyc2','skyc3'};

%-------------------------------------------------------------------------%
% Constants

% Conversions
mph_to_mps = 0.44704;       % miles-per-hour to meters-per-second
mma_to_cms = 10/(24*60*60); % 1-hour accumulation to centimeters-per-second

% Station info
wthr_info = {'BID','BLOCK ISLAND (AWOS)', 41.17,   -71.58;    ...   % RI
             'UUU','NEWORT',              41.53244,-71.28154; ...   % RI
             'OQU','N.KINGSTON/QUONSET',  41.59714,-71.41214; ...   % RI
             'SVZ','Pawtucket (AWOS)',    41.92076,-71.49138; ...   % RI
             'PVD','PROVIDENCE/GREEN',    41.7219, -71.4325;  ...   % RI
             'WST','WESTERLY',            41.3497, -71.7989;  ...   % RI
             'GON','GROTON/NEW LONDON',   41.33,   -72.05;    ...   % CT
             'HVN','NEW HAVEN/TWEED',     41.26375,-72.88681; ...   % CT
             'TAN','TAUNTON MUNI AIRPORT',41.87556,-71.02111; ...   % MA
             'EWB','NEW BEDFORD MUNI',    41.67639,-70.95833; ...   % MA
             'MVY','MARTHAS VINEYARD',    41.39306,-70.615;   ...   % MA/IS
             'ACK','NANTUCKET MEMORIAL',  41.25311,-70.06031; ...   % MA/IS
             'CQX','CHATHAM MUNI ARPT',   41.6875, -69.99333; ...   % MA
             'HYA','HYANNIS/POLANDO FLD'  41.67,   -70.28;    ...   % MA
             'MTP','MONTAUK AIRPORT',     41.07306,-71.92333; ...   % NY/LI
             'HTO','EAST HAMPTON',        40.96,   -72.25;    ...   % NY/LI
             'FOK','WESTHAMPTON BEACH',   40.84365,-72.63179; ...   % NY/LI
            };
wthr_to_get = wthr_info(:,1);

%=========================================================================%
% Plot of station locations

% la = shaperead('landareas.shp');
% plot([la.X],[la.Y],'-k');
% hold on;
% for i=1:size(wthr_info,1)
%     plot(wthr_info{i,4},wthr_info{i,3},'sr');
%     text(wthr_info{i,4},wthr_info{i,3}+0.01,wthr_info{i,1});
% end
% hold off;
% xlim([-73 -69]); ylim([40.5 42.2]);

%=========================================================================%

% Loop through stations
n_wthr = numel(wthr_to_get);
n_var  = numel(var_to_get);
n_year = numel(year_to_get);
for iy=1:n_year
for iw=1:n_wthr
clc; disp(['On year ' num2str(year_to_get(iy)) ', Station ' wthr_to_get{iw} '...']);

    % Create save directory if it doesn't yet exist
    wthr_dir = [save_dir '/' wthr_to_get{iw}];
    if( exist(wthr_dir,'dir')~=7 )
        mkdir(wthr_dir);
    end
    
    % Continue only if data file does not yet exist
    save_file = [wthr_dir '/' wthr_to_get{iw} '_' num2str(year_to_get(iy)) '.mat'];
    if(exist(save_file,'file')~=2)
    
        % Construct data URL
        wthr_url = ['https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=' wthr_to_get{iw}];
        for iv=1:numel(var_to_get); wthr_url = [wthr_url '&data=' var_to_get{iv}]; end
        if(any(strcmp(var_to_get,'skyc1')==1)); wthr_url = [wthr_url '&data=metar'];    end
        wthr_url = [wthr_url '&year1=' num2str(year_to_get(iy)) '&month1=1&day1=1'   ...
                             '&year2=' num2str(year_to_get(iy)) '&month2=12&day2=31' ...
                    '&latlon=no&elev=no&missing=empty&trace=0.0001&direct=no&report_type=1&report_type=2'];         

        % Load data within specified time span
        wthr_raw = webread(wthr_url,weboptions('ContentType','text','Timeout',inf));

        % Save as .csv file and read in as a table
        fid = fopen('wthr_raw.csv','wt');
        fprintf(fid,wthr_raw);
        fclose(fid);
        wthr_table = readtable('wthr_raw.csv');
        delete('wthr_raw.csv');

        % If now empty, there's no data available
        if(isempty(wthr_table))
            warning(['No data exists for station ' wthr_to_get{iw} ' for the specified time frame, skipping...']);
        else

            % Get available variables
            wthr_vars = wthr_table.Properties.VariableNames;

            % Get time variable
            wthr_time = datenum(wthr_table.valid);
            nt        = numel(wthr_time);

            % Variables to save to file
            save_str = {};

            % Gather variables
            for iv=1:numel(wthr_vars)
            switch wthr_vars{iv}
                case 'tmpc'
                if(~all(isnan(wthr_table.tmpc)))
                    Tair      = wthr_table.tmpc;    
                    %tair_time = wthr_time;
                    %Tair      = nanfill(Tair);
                    save_str{end+1} = 'Tair';   %:end+2) = {'Tair','tair_time'};
                end
                case 'mslp'   
                if(~all(isnan(wthr_table.mslp)))
                    Pair      = wthr_table.mslp;    
                    %pair_time = wthr_time;
                    %Pair      = nanfill(Pair);
                    save_str{end+1} = 'Pair';   %(end+1:end+2) = {'Pair','pair_time'};
                end
                case 'relh'
                if(~all(isnan(wthr_table.relh)))
                    Qair      = wthr_table.relh;    
                    %qair_time = wthr_time;
                    %Qair      = nanfill(Qair);
                    save_str{end+1} = 'Qair';   %(end+1:end+2) = {'Qair','qair_time'};
                end
                case 'p01m'
                if(~all(isnan(wthr_table.p01m)))
                    rain              = wthr_table.p01m .* mma_to_cms;    
                    %rain_time         = wthr_time;
                    %rain(isnan(rain)) = 0;
                    save_str{end+1} = 'rain'; %(end+1:end+2) = {'rain','rain_time'};
                end
                case 'drct'
                if(~all(isnan(wthr_table.drct)))
                    Uwind     = wthr_table.sped .* sin(wthr_table.drct) .* mph_to_mps;    
                    Vwind     = wthr_table.sped .* cos(wthr_table.drct) .* mph_to_mps;   
                    iwind     = find(~isnan(Uwind) & ~isnan(Vwind));
                    if(~isempty(iwind))
                        Zwind = Uwind + sqrt(-1).*Vwind;
                        %Zwind = nanfill(Zwind);
                        Uwind(iwind) = real(Zwind(iwind));
                        Vwind(iwind) = imag(Zwind(iwind));
                    end
                    %wind_time = wthr_time;
                    save_str(end+1:end+2) = {'Uwind','Vwind'}; %,'wind_time'};
                    clear iwind;
                end
                case 'skyc1'
                    if(~iscell(wthr_table.skyc1));  tmp1 = cell(numel(wthr_table.skyc1),1); else;   tmp1=wthr_table.skyc1;  end
                    if(~iscell(wthr_table.skyc2));  tmp2 = cell(numel(wthr_table.skyc2),1); else;   tmp2=wthr_table.skyc2;  end
                    if(~iscell(wthr_table.skyc3));  tmp3 = cell(numel(wthr_table.skyc3),1); else;   tmp3=wthr_table.skyc3;  end
                    Cstr = [tmp1, tmp2, tmp3];
                    Mstr = wthr_table.metar;
                    Cfra = NaN(nt,1);
                    for it=1:nt
                    if(~all(isempty([Cstr{it,:}])))

                       % Available cloud values
                       tmp = NaN(3,1);
                       for i=1:3
                       if(~isempty(Cstr{it,i}))    
                       switch Cstr{it,i}
                           case 'FEW';  tmp(i) = 1.5/8;
                           case 'SCT';  tmp(i) = 3.5/8;
                           case 'BKN';  tmp(i) = 5.5/8;
                           case 'OVC';  tmp(i) = 1;
                       end
                       end
                       end
                       
                       % Take max value
                       Cfra(it) = nanmax(tmp);
                       clear tmp;

                    end
                    end
                    clear it Cstr Mstr;

                    % Save
                    %cloud_time = wthr_time;
                    %Cfra       = nanfill(Cfra);
                    save_str{end+1} = 'Cfra';   %(end+1:end+2) = {'Cfra','cloud_time'};
            end    
            end
            clear iv;
            
            % Save data
            time = wthr_time;
            eval_str = 'save(save_file,';
            for i=1:numel(save_str);    eval_str = [eval_str, '''' save_str{i} ''',']; end
            eval_str = [eval_str '''time'');'];
            eval(eval_str);
            %save_str{end+1} = 'time';
            %save(save_file,save_str);

            % Clean-up
            clear nt Tair Pair Qair rain Uwind Vwind Cfra save_str *_time time;

        end
        clear wthr_url wthr_raw wthr_table;
    
    end
    clear wthr_dir save_file;
    
end
end
clear iw iy;

%=========================================================================%