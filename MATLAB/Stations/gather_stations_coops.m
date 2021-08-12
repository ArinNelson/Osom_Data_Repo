%=========================================================================%
% Gather land/dock station data from CO-OPS
% Arin Nelson
% on 07/23/2021
% 
% Last edited 08/01/2021
%=========================================================================%
clear; clc; clear mex; addpath('../Utilities');

% Options
get_year = 1990:2020;
save_dir = 'D:/OSOM_Data_Repo/Stations/COOPS/';
%sttn_id  = {'8452944','8454049','8452660','8447386','8447386','8447386',...
%            '8447386','8447930','8447435','8449130','8461490','8465705',...
%            '8467150','8510560'};
sttn_id = {'8454000','8452944','8454049','8452660','8447386'};

% Constants
var_to_get = {'water_level','water_temperature','conductivity','salinity',...
              'air_temperature','wind','air_pressure'}; % Variables to look for
n_var = numel(var_to_get);          

%=========================================================================%
for is=1:numel(sttn_id)

    % Station directory
    sttn_savedir = [save_dir sttn_id{is}];
    if(exist(sttn_savedir,'dir')~=7)
        mkdir(sttn_savedir);
    end

    % Station info file
    sttn_savefile = [sttn_savedir '\' sttn_id{is} '_info.mat'];
    if(exist(sttn_savefile,'file')~=2)

        % Get station information
        sttn_url        = ['https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' sttn_id{is} '.xml'];
        xml_raw         = xmlread_web(sttn_url);
        sttn_raw        = xml_raw(1).Stations(1).Station(1);
        name            = sttn_raw.name.Text;
        lat             = str2double( sttn_raw.lat.Text );
        lon             = str2double( sttn_raw.lng.Text );
        affiliations    = sttn_raw.affiliations.Text;
        is_tidal        = sttn_raw.tidal.Text;
        state           = sttn_raw.state.Text;
        timezone_offset = str2double( sttn_raw.timezonecorr.Text );
        clear sttn_url xml_raw sttn_raw;

        % Get station tides
        sttn_url = ['https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/' sttn_id{is} '/harcon.xml?units=metric'];xml_raw  = xmlread_web(sttn_url);
        sttn_raw = xml_raw(1).HarmonicConstituents(1).HarmonicConstituent;
        tide     = struct;
        for it=1:numel(sttn_raw)
            tide(it).name      = sttn_raw{it}.name.Text;
            tide(it).amplitude = str2double(sttn_raw{it}.amplitude.Text);
            tide(it).phase_GMT = str2double(sttn_raw{it}.phase_GMT.Text); 
        end
        clear sttn_url xml_raw sttn_raw it;

        % Save
        save(sttn_savefile,'name','lat','lon','affiliations','is_tidal','state','timezone_offset','tide');
        clear name lat lon affiliations is_tidal state timezone_offset tide;

    end

    %---------------------------------------------------------------------%

    % Loop through years & months
    for iy=1:numel(get_year)
    for im=1:12
    clc; disp(['Gathering data for station ' sttn_id{is} ', year/month ' num2str(get_year(iy)) '/' sprintf('%0.2d',im) '...']);

        % Base url string
        url_base = ['https://api.tidesandcurrents.noaa.gov/api/prod/datagetter' ...
                    '?begin_date='  num2str(get_year(iy)) sprintf('%0.2d',im) '01' ...
                    '&end_date=' num2str(get_year(iy)) sprintf('%0.2d',im) num2str(eomday(get_year(iy),im)) ...
                    '&station=' sttn_id{is} ...
                    '&datum=STND&time_zone=gmt&units=metric&format=csv&product='];

        % Loop through variables
        for iv=1:n_var

            % Save directory
            var_savedir = [sttn_savedir '\' var_to_get{iv}];
            if(exist(var_savedir,'dir')~=7);   mkdir(var_savedir); end

            % Save file
            var_savefile = [var_savedir '\' sttn_id{is} '_' var_to_get{iv} '_' num2str(get_year(iy)) '_' sprintf('%0.2d',im) '.mat'];
            if(exist(var_savefile,'file')~=2)

                % Try to read the url
                try txt_raw = webread([url_base var_to_get{iv}]);
                catch err
                  txt_raw = webread([url_base var_to_get{iv}],weboptions('ContentType','text'));
                  txt_splt = strsplit(txt_raw,{'\n'});
                  txt_splt = txt_splt(2:end-1);
                  txt_raw  = cell(0,2);
                  for i=1:numel(txt_splt)
                    txt_line = strsplit(txt_splt{i},',');
                    if(~isempty(txt_line{2}))
                      txt_raw{end+1,1} = txt_line(1);
                      txt_raw{end,2}   = str2double(txt_line{2});
                    end
                  end
                  clear txt_splt txt_line i;
                end
                if(~isempty(txt_raw))

                    % Parse the data
                    time  = datenum([txt_raw{:,1}]);
                    value = txt_raw{:,2};

                    % Save the data to file
                    save(var_savefile,'time','value');

                    % Clean-up
                    clear time value;

                end
                clear txt_raw;

            end
            clear var_savefile;

        end
        clear iv;

        % Clean-up
        clear url_base;

    end
    end
    clear iy im;

end
%=========================================================================%