% Download data for USGS streamgages and compute the daily statistics
ggl_id     = '1Uc7mBJH-hggEDFPP6eWMh6tITHHFbDG2jtcLnI0ubdo';
ggl_sheet  = 'Fixed_Obs';
ggl_url    = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',ggl_id,ggl_sheet);
tbl_obs    = webread(ggl_url);
i_nws      = find( strcmp( tbl_obs.Type, 'Weather Station' )==1 & strcmp(tbl_obs.Use,'T')==1 );
nws_id     = tbl_obs.ID(i_nws);
start_date = '1990-01-01';
save_dir   = 'F:/OSOM_Data_Repo/Stations/NWS/';

% Gather data
for i=1:numel(nws_id)
clc; disp(['On weather station # ' num2str(i) ' of ' num2str(numel(nws_id)) '...']);
    download_nws_history(nws_id{i},[save_dir 'nws_' nws_id{i} '_raw.nc'],start_date);
end
clear i;

