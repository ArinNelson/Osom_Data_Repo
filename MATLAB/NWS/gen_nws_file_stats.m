function gen_nws_file_stats(nws_id,save_file,varargin)

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
    dimid = [];
    dimid(end+1) = netcdf.defDim(ncid,'time',nc_unlimited);
    dimid(end+1) = netcdf.defDim(ncid,'percentile',7);
    
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
    varid(end+1) = netcdf.defVar(ncid,'time','NC_DOUBLE',dimid(1));
        netcdf.putAtt(ncid,varid(end),'units','matlab datenum');
        netcdf.putAtt(ncid,varid(end),'timezone','UTC');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    
    % Dynamic variables
    str_varname = {'tmpc','relh','drct','sped','mslp','p01m','cfra'};
    str_varlong = {'Air Temperature','Relative Humidity','Wind Direction','Wind Speed','Sea Level Pressure','Precipitation Rate','Cloud Fraction'};
    str_varunit = {'degrees_Celsius','percent','degrees_from_true_north','miles_per_hour','millibar','inches_per_day','percent'};
    for iv=1:numel(str_varname)
        varid(end+1) = netcdf.defVar(ncid,[str_varname{iv} '_mean'],'NC_DOUBLE',dimid(1));
            netcdf.putAtt(ncid,varid(end),'long_name',['Daily Mean ' str_varlong{iv}]);
            netcdf.putAtt(ncid,varid(end),'units',str_varunit{iv});
            netcdf.defVarDeflate(ncid,varid(end),true,true,2);
            netcdf.defVarFill(ncid,varid(end),false,-999);
        varid(end+1) = netcdf.defVar(ncid,[str_varname{iv} '_std'],'NC_DOUBLE',dimid(1));
            netcdf.putAtt(ncid,varid(end),'long_name',['Daily Standard Deviation of ' str_varlong{iv}]);
            netcdf.putAtt(ncid,varid(end),'units',str_varunit{iv});
            netcdf.defVarDeflate(ncid,varid(end),true,true,2);
            netcdf.defVarFill(ncid,varid(end),false,-999);
        varid(end+1) = netcdf.defVar(ncid,[str_varname{iv} '_num'],'NC_DOUBLE',dimid(1));
            netcdf.putAtt(ncid,varid(end),'long_name',['Daily Number of ' str_varlong{iv} ' Observations']);
            netcdf.putAtt(ncid,varid(end),'units','count');
            netcdf.defVarDeflate(ncid,varid(end),true,true,2);
            netcdf.defVarFill(ncid,varid(end),false,-999);
        varid(end+1) = netcdf.defVar(ncid,[str_varname{iv} '_cvg'],'NC_DOUBLE',dimid(1));
            netcdf.putAtt(ncid,varid(end),'long_name',['Percent Coverage of Daily ' str_varlong{iv} ' Observations Captured']);
            netcdf.putAtt(ncid,varid(end),'units','percent');
            netcdf.defVarDeflate(ncid,varid(end),true,true,2);
            netcdf.defVarFill(ncid,varid(end),false,-999);
        varid(end+1) = netcdf.defVar(ncid,[str_varname{iv} '_prctile'],'NC_DOUBLE',dimid([1 2]));
            netcdf.putAtt(ncid,varid(end),'long_name',['Daily Percentiles of ' str_varlong{iv}]);
            netcdf.putAtt(ncid,varid(end),'units',str_varunit{iv});
            netcdf.defVarDeflate(ncid,varid(end),true,true,2);
            netcdf.defVarFill(ncid,varid(end),false,-999);
    end
         
    % End definitions
    netcdf.endDef(ncid);
    netcdf.close(ncid);
    
    % Write station metadata
    ncwrite(save_file,'longitude',lon );
    ncwrite(save_file,'latitude', lat );
    ncwrite(save_file,'elevation',elev);
        
end