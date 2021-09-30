%function download_riddc_landsat
%=========================================================================%
% (header)
%=========================================================================%

% RIDDC data url's
landsat_grid_url = 'https://pricaimcit.services.brown.edu/erddap/griddap/landsat_sst_fc65_afa3_6231';
landsat_data_url = 'https://pricaimcit.services.brown.edu/erddap/griddap/landsat_sst_572a_70ac_7b1d';

% Notes from netcdf attributes
% 
% Attribute Accuracy Report: Satellite-derived orthorectified brightness 
% temperature was measured within 0.1 degrees C for Landsat 8, 0.6 degrees 
% C for Landsat 7, and 0.5 degrees C for Landsat 5. Satellite measurements 
% were compared to in situ (buoy) surface temperatures from 2003 to 2015, 
% and the mean bias between the RI DEM buoy temperatures and the satellite 
% temperatures at the pixel of the buoys was added to all scenes by 
% satellite. The standard deviation between the buoys and the satellie 
% after adding the bias is 1.9 for Landsat 5, 1.9 for Landsat 7 and 1.3 for 
% Landsat 8. The standard deviation is considered the uncertainty of the 
% satellite measurements. See https://www.usgs.gov/land-resources/nli/
% landsat/landsat-surface-reflectance?qt-science_support_page_related_con=0
% #qt-science_support_page_related_con for more information.
% 
% Converted from Landsat 5, 7, and 8 Surface Refelectance geotiff products 
% to netCDF. The units were changed from K to degrees C and the average 
% bias determined through RI DEM buoy comparison in Narragansett Bay 
% (2003-2015) to each satellite was added to all scenes from the 
% corresponding satellite. For Landsat 5, 3.36 degrees C was added to all 
% scenes, 3.34 for Landsat 7, and 1.92 for Landsat 8. The errors were 
% determined by calculating the standard deviation of the difference 
% between Landsat and buoy temperature at buoy locations for each 
% satellite. The scenes were also cloud masked, land masked, and stripes in 
% Landsat 7 imagery (due to sensor failure) were masked as well. Data 
% outside of the Narragansett Bay region and for all cloud cover is 
% included, though a buoy comparison was only conducted within Narragansett 
% Bay for Landsat scenes with less than 50% cloud cover. As a result, data 
% uncertainties are unknown outside of the Narragansett Bay region and for 
% scenes with greater cloud cover.

save_file = 'F:\OSOM_Data_Repo\LandSat\landsat_riddc.nc';

%-------------------------------------------------------------------------%

% File info's
grid_info = ncinfo(landsat_grid_url);
data_info = ncinfo(landsat_data_url);

% Load grid
x   = ncread(landsat_grid_url,'X');
y   = ncread(landsat_grid_url,'Y'); 
lon = ncread(landsat_grid_url,'Longitude');
lat = ncread(landsat_grid_url,'Latitude');

% 1D data
time     = ncread(landsat_data_url,'time');          % seconds since 1970-01-01T00:00:00Z
cfra_all = squeeze(ncread(landsat_data_url,'Clouds',   [1 1 1],[inf 1 1]));
cfra_nb  = squeeze(ncread(landsat_data_url,'NG_Clouds',[1 1 1],[inf 1 1]));
satnum   = squeeze(ncread(landsat_data_url,'satellite',[1 1 1],[inf 1 1]));

%-------------------------------------------------------------------------%

% NetCDF variables
nc_global    = netcdf.getConstant('NC_GLOBAL');
nc_unlimited = netcdf.getConstant('NC_UNLIMITED');

% Initialize the NetCDF file
ncid = netcdf.create(save_file,'NETCDF4');

% Save the global attributes
for i=1:numel(data_info.Attributes)
    netcdf.putAtt(ncid,nc_global,data_info.Attributes(i).Name,data_info.Attributes(i).Value);
end

% Define the dimensions
dimid        = [];
dimid(end+1) = netcdf.defDim(ncid,'x',numel(x));
dimid(end+1) = netcdf.defDim(ncid,'y',numel(y));
dimid(end+1) = netcdf.defDim(ncid,'t',nc_unlimited);    

% Define the variables
varid        = [];
varid(end+1) = netcdf.defVar(ncid,'longitude','NC_DOUBLE',dimid([1 2]));
        netcdf.putAtt(ncid,varid(end),'units','degrees_east');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
varid(end+1) = netcdf.defVar(ncid,'latitude','NC_DOUBLE',dimid([1 2]));
        netcdf.putAtt(ncid,varid(end),'units','degrees_north');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);   
varid(end+1) = netcdf.defVar(ncid,'time','NC_DOUBLE',dimid(3));
        netcdf.putAtt(ncid,varid(end),'units','matlab_datenum');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);  
varid(end+1) = netcdf.defVar(ncid,'cfra_all','NC_DOUBLE',dimid(3));
        netcdf.putAtt(ncid,varid(end),'long_name','cloud_fraction_entire_domain');
        netcdf.putAtt(ncid,varid(end),'units','percent');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);     
varid(end+1) = netcdf.defVar(ncid,'cfra_nb','NC_DOUBLE',dimid(3));
        netcdf.putAtt(ncid,varid(end),'long_name','cloud_fraction_narragansett_bay_only');
        netcdf.putAtt(ncid,varid(end),'units','percent');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2); 
varid(end+1) = netcdf.defVar(ncid,'sat_num','NC_INT',dimid(3));
        netcdf.putAtt(ncid,varid(end),'long_name','landsat_satellite_number');
        netcdf.putAtt(ncid,varid(end),'units','integer');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);  
varid(end+1) = netcdf.defVar(ncid,'temp','NC_DOUBLE',dimid([1 2 3]));
        netcdf.putAtt(ncid,varid(end),'long_name','surface_temperature');
        netcdf.putAtt(ncid,varid(end),'Accuracy_of_Satellite','within 0.1 degrees C for Landsat 8, 0.6 degrees C for Landsat 7, and 0.5 degrees C for Landsat 5. Uncertainty measurment of 1.9 degrees C derived from Landsat 5 comaprison with RI DEM buoys (2003-2011).');
        netcdf.putAtt(ncid,varid(end),'Measurement_Uncertainty_derived_from_surface_buoy_comparison_degC','1.3');
        netcdf.putAtt(ncid,varid(end),'units','degrees_celsius');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);          
    
% End definitions
netcdf.endDef(ncid);
netcdf.close(ncid);
    
% Save 1D variables
ncwrite(save_file,'longitude',lon);
ncwrite(save_file,'latitude',lat);
ncwrite(save_file,'time',(time./(60*60*24))+datenum('1970-01-01'));
ncwrite(save_file,'cfra_all',cfra_all);
ncwrite(save_file,'cfra_nb',cfra_nb);
ncwrite(save_file,'sat_num',satnum);

%-------------------------------------------------------------------------%

% Download 2D data
for it=122:numel(time)
clc; disp(['on timestep ' num2str(it) '/' num2str(numel(time)) '...']);

    % Read
	temp = squeeze( ncread(landsat_data_url,'Temp',[it 1 1],[1 inf inf]) );
    
    % Save
    ncwrite(save_file,'temp',temp,[1 1 it]);
    
end
clear it temp;

%-------------------------------------------------------------------------%

% Make plots of the data
x  = ncread(save_file,'longitude');
y  = ncread(save_file,'latitude');
t  = ncread(save_file,'time');

% Keep data for 2018+
it = find(year(t)<2018,1,'last') : numel(t);
nt = numel(it);

% Loop
for i=1:nt

    % Get data
    sst = ncread(save_file,'temp',    [1 1 it(i)],[inf inf 1]);
    ca  = ncread(save_file,'cfra_all',it(i),1);
    cnb = ncread(save_file,'cfra_nb', it(i),1);
    
    % Generate plot
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    axx = axes('units','normalized','position',[0.025 0.05 0.95 0.90],'box','on');
    
    % Plot
    surf(x,y,sst,'edgecolor','none'); view(2); caxis([0 25]); colormap(jet);
    xlabel('longitude'); ylabel('latitude');
    cb = colorbar; set(get(cb,'ylabel'),'string','SST (^oC)');
    title(['LandSat SST, ' datestr(t(it(i))) ', ' sprintf('%0.3f',cnb) '/' sprintf('%0.3f',ca)]);
    set(gca,'xlim',[-72.87 -69.82],'ylim',[40.66 42.2],'color','k');
    
    % Save plot
    print(gcf,['Plots/LandSatSST_' sprintf('%0.3d',it(i)) '.png'],'-dpng');
    
    % Clean-up
    close all; clear sst ca cnb fig axx cb;
    
end
clear i;

%-------------------------------------------------------------------------%