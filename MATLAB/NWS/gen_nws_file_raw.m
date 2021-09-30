function gen_nws_file_raw(nws_id,save_file,varargin)

%-------------------------------------------------------------------------%
% GATHER SITE METADATA & VARIABLES IN INVENTORY

    % to do...


%-------------------------------------------------------------------------%
% CREATE THE NETCDF FILE      
    
    % Init the NetCDF
    if(exist(save_file,'file')==2); delete(save_file);  end
    ncid = netcdf.create(save_file,'NETCDF4');
    
    % NetCDF constants
    nc_global    = netcdf.getConstant('NC_GLOBAL');
    nc_unlimited = netcdf.getConstant('NC_UNLIMITED');
    
    % Define the global attributes
    netcdf.putAtt(ncid,nc_global,'title',  ['NWS weather station historical data file (raw), ID # ' nws_id]);
    netcdf.putAtt(ncid,nc_global,'history',['Created by gen_nws_file_raw.m on ' datestr(now)]);
    
    % Titles are saved internally for now...
    station_name = '';
    switch nws_id
        case {'BID','KBID'};    station_name = 'BLOCK ISLAND (AWOS)';   station_network = 'RI_ASOS';    lat = 41.17;        lon = -71.58;       elev = 33; 
        case {'UUU','KUUU'};    station_name = 'NEWPORT';               station_network = 'RI_ASOS';    lat = 41.53244;     lon = -71.28154;    elev = 29;   
        case {'OQU','KOQU'};    station_name = 'N. KINGSTON/QUONSET';   station_network = 'RI_ASOS';    lat = 41.59714;     lon = -71.41214;    elev = 6;
        case {'SFZ','KSFZ'};    station_name = 'PAWTUCKET (AWOS)';      station_network = 'RI_ASOS';    lat = 41.92076;     lon = -71.49138;    elev = 134;
        case {'PVD','KPVD'};    station_name = 'PROVIDENCE/GREEN';      station_network = 'RI_ASOS';    lat = 41.7219;      lon = -71.4325;     elev = 19;
        case {'WST','KWST'};    station_name = 'Westerly';              station_network = 'RI_ASOS';    lat = 41.3497;      lon = -71.7989;     elev = 24;
        case {'EWB','KEWB'};    station_name = 'NEW BEDFORD MUNI';      station_network = 'MA_ASOS';    lat = 41.67639;     lon = -70.95833;    elev = 24;
        case {'TAN','KTAN'};    station_name = 'TAUNTON MUNI AIRPORT';  station_network = 'MA_ASOS';    lat = 41.87556;     lon = -71.02111;    elev = 13;
        case {'GON','KGON'};    station_name = 'GROTON/NEW LONDON';     station_network = 'CT_ASOS';    lat = 41.33;        lon = -72.05;       elev = 3;
        case {'HFD','KHFD'};    station_name = 'HARTFORD/BRAINARD';     station_network = 'CT_ASOS';    lat = 41.73672;     lon = -72.64944;    elev = 6;
    end
    netcdf.putAtt(ncid,nc_global,'station_id',nws_id);
    netcdf.putAtt(ncid,nc_global,'station_name',station_name);
    netcdf.putAtt(ncid,nc_global,'station_network',station_network);
    
    % Define the time dimension
    dimid = netcdf.defDim(ncid,'time',nc_unlimited);
    
    % Define variables
    varid = [];
    
    % Location data
    varid(end+1) = netcdf.defVar(ncid,'longitude','NC_DOUBLE',[]);
        netcdf.putAtt(ncid,varid(end),'units','degrees_east');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    varid(end+1) = netcdf.defVar(ncid,'latitude','NC_DOUBLE',[]);
        netcdf.putAtt(ncid,varid(end),'units','degrees_north');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    varid(end+1) = netcdf.defVar(ncid,'elevation','NC_DOUBLE',[]);
        netcdf.putAtt(ncid,varid(end),'units','meters_above_msl');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
        
    % Time variable
    varid(end+1) = netcdf.defVar(ncid,'time','NC_DOUBLE',dimid);
        netcdf.putAtt(ncid,varid(end),'units','matlab datenum');
        netcdf.putAtt(ncid,varid(end),'timezone','UTC');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    
    % Dynamic variables
    varid(end+1) = netcdf.defVar(ncid,'tmpc','NC_DOUBLE',dimid);
        netcdf.putAtt(ncid,varid(end),'long_name','Air Temperature');
        netcdf.putAtt(ncid,varid(end),'units','degrees_Celsius');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    varid(end+1) = netcdf.defVar(ncid,'relh','NC_DOUBLE',dimid);
        netcdf.putAtt(ncid,varid(end),'long_name','Relative Humidity');
        netcdf.putAtt(ncid,varid(end),'units','percent');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    varid(end+1) = netcdf.defVar(ncid,'drct','NC_DOUBLE',dimid);
        netcdf.putAtt(ncid,varid(end),'long_name','Wind Direction');
        netcdf.putAtt(ncid,varid(end),'units','degrees_from_true_north');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    varid(end+1) = netcdf.defVar(ncid,'sped','NC_DOUBLE',dimid);
        netcdf.putAtt(ncid,varid(end),'long_name','Wind Speed');
        netcdf.putAtt(ncid,varid(end),'units','miles_per_hour');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    varid(end+1) = netcdf.defVar(ncid,'mslp','NC_DOUBLE',dimid);
        netcdf.putAtt(ncid,varid(end),'long_name','Sea Level Pressure');
        netcdf.putAtt(ncid,varid(end),'units','millibar');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    varid(end+1) = netcdf.defVar(ncid,'p01m','NC_DOUBLE',dimid);
        netcdf.putAtt(ncid,varid(end),'long_name','One Hour Precipitation');
        netcdf.putAtt(ncid,varid(end),'units','inches_since_last_reset');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    varid(end+1) = netcdf.defVar(ncid,'cfra','NC_DOUBLE',dimid);
        netcdf.putAtt(ncid,varid(end),'long_name','Cloud Fraction');
        netcdf.putAtt(ncid,varid(end),'units','percent');
        netcdf.putAtt(ncid,varid(end),'method','max(skyc1,skyc2,skyc3)');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
         
    % End definitions
    netcdf.endDef(ncid);
    netcdf.close(ncid);
    
    % Write station metadata
    ncwrite(save_file,'longitude',lon );
    ncwrite(save_file,'latitude', lat );
    ncwrite(save_file,'elevation',elev);
        
end