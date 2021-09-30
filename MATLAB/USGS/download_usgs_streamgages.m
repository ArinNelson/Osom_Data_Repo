% Download data for USGS streamgages and compute the daily statistics
ggl_id     = '1Uc7mBJH-hggEDFPP6eWMh6tITHHFbDG2jtcLnI0ubdo';
ggl_sheet  = 'Fixed_Obs';
ggl_url    = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',ggl_id,ggl_sheet);
tbl_obs    = webread(ggl_url);
i_rvr      = find( strcmp( tbl_obs.Type, 'River Gauge' )==1 & strcmp(tbl_obs.Use,'T')==1 );
usgs_id    = tbl_obs.ID(i_rvr);
start_date = tbl_obs.Start(i_rvr);
save_dir   = 'F:/OSOM_Data_Repo/USGS/RiverGauges/';

% Gather data
for i=16:numel(usgs_id)
clc; disp(['On river # ' num2str(i) ' of ' num2str(numel(usgs_id)) '...']);
    save_file = [save_dir 'usgs_' usgs_id{i} '_raw.nc'];
    sd        = start_date{i}([1:4 6:7 9:10]);
    download_usgs_streamgage(usgs_id{i},save_file,sd);
end
clear i save_file sd;

% Compute daily statistics
parfor i=16:numel(usgs_id)
    process_usgs_streamgage([save_dir 'usgs_' usgs_id{i} '_raw.nc'],[save_dir 'usgs_' usgs_id{i} '_daily.nc']);
end
clear i;