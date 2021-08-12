function status = nc_gen_frc_roms(opts,nx,ny,var_to_incl)

    % netCDF constants
    nc_global    = netcdf.getConstant('NC_GLOBAL');
    nc_unlimited = netcdf.getConstant('NC_UNLIMITED');

    % NetCDF creation mode options
    ncopts = netcdf.getConstant('NETCDF4');

    % Create the NetCDF   
    if(exist(opts.FileName,'file')==2); delete(opts.FileName); end
    ncid = netcdf.create(opts.FileName,ncopts);
  
    % Global attributes
    netcdf.putAtt(ncid,nc_global,'type','FORCING file');
    netcdf.putAtt(ncid,nc_global,'title',opts.Title);
    netcdf.putAtt(ncid,nc_global,'history',['generated by nc_gen_frc_roms.m on ' datestr(now)]);
  
    % Grid dimensions
    dimid = [];
    dimid(end+1) = netcdf.defDim(ncid,'lon',nx);
    dimid(end+1) = netcdf.defDim(ncid,'lat',ny);
    for i=1:numel(var_to_incl)
    switch var_to_incl{i}
        case 'winds';   dimid(end+1) = netcdf.defDim(ncid,'wind_time',nc_unlimited);
        case 'Pair';    dimid(end+1) = netcdf.defDim(ncid,'pair_time',nc_unlimited);
        case 'Tair';    dimid(end+1) = netcdf.defDim(ncid,'tair_time',nc_unlimited);  
        case 'Qair';    dimid(end+1) = netcdf.defDim(ncid,'qair_time',nc_unlimited);
        case 'Cfra';    dimid(end+1) = netcdf.defDim(ncid,'cloud_time',nc_unlimited);
        case 'rain';    dimid(end+1) = netcdf.defDim(ncid,'rain_time',nc_unlimited);
        case {'lwrad','lwrad_down'};	dimid(end+1) = netcdf.defDim(ncid,'lrf_time',nc_unlimited);
        case {'swrad','swrad_down'};	dimid(end+1) = netcdf.defDim(ncid,'srf_time',nc_unlimited);    
    end
    end
  
    % Grid variables
    varid = [];
    varid(end+1) = netcdf.defVar(ncid,'lon','double',dimid(1));
        netcdf.putAtt(ncid,varid(end),'long_name','longitude');
        netcdf.putAtt(ncid,varid(end),'units','degrees east');
    varid(end+1) = netcdf.defVar(ncid,'lat','double',dimid(1));
        netcdf.putAtt(ncid,varid(end),'long_name','latitude');
        netcdf.putAtt(ncid,varid(end),'units','degrees north');
  
    % Data Variables
    for i=1:numel(var_to_incl)
    switch var_to_incl{i}
    
        %-----------------------------------------------------------------%  
        case 'winds'
          varid(end+1) = netcdf.defVar(ncid,'wind_time','double',dimid(i+2));
            netcdf.putAtt(ncid,varid(end),'long_name','surface wind time');
            netcdf.putAtt(ncid,varid(end),'units','days since initialization');
          varid(end+1) = netcdf.defVar(ncid,'Uwind','double',[dimid(1) dimid(2) dimid(i+2)]);
            netcdf.putAtt(ncid,varid(end),'long_name','surface u-wind component');
            netcdf.putAtt(ncid,varid(end),'units','meters second-1');
            netcdf.putAtt(ncid,varid(end),'time','wind_time');
            netcdf.putAtt(ncid,varid(end),'coordinates','lon lat wind_time');
          varid(end+1) = netcdf.defVar(ncid,'Vwind','double',[dimid(1) dimid(2) dimid(i+2)]);
            netcdf.putAtt(ncid,varid(end),'long_name','surface v-wind component');
            netcdf.putAtt(ncid,varid(end),'units','meters second-1');
            netcdf.putAtt(ncid,varid(end),'time','wind_time');
            netcdf.putAtt(ncid,varid(end),'coordinates','lon lat wind_time');

        %-----------------------------------------------------------------%  
        case 'Pair'
          varid(end+1) = netcdf.defVar(ncid,'pair_time','double',dimid(i+2));
            netcdf.putAtt(ncid,varid(end),'long_name','surface air pressure time');
            netcdf.putAtt(ncid,varid(end),'units','days since initialization');
          varid(end+1) = netcdf.defVar(ncid,'Pair','double',[dimid(1) dimid(2) dimid(i+2)]);
            netcdf.putAtt(ncid,varid(end),'long_name','surface air pressure');
            netcdf.putAtt(ncid,varid(end),'units','millibars');
            netcdf.putAtt(ncid,varid(end),'time','pair_time');
            netcdf.putAtt(ncid,varid(end),'coordinates','lon lat pair_time');

        %-----------------------------------------------------------------%  
        case 'Tair'
          varid(end+1) = netcdf.defVar(ncid,'tair_time','double',dimid(i+2));
            netcdf.putAtt(ncid,varid(end),'long_name','surface air temperature time');
            netcdf.putAtt(ncid,varid(end),'units','days since initialization');
          varid(end+1) = netcdf.defVar(ncid,'Tair','double',[dimid(1) dimid(2) dimid(i+2)]);
            netcdf.putAtt(ncid,varid(end),'long_name','surface air temperature');
            netcdf.putAtt(ncid,varid(end),'units','Celsius');
            netcdf.putAtt(ncid,varid(end),'time','tair_time');
            netcdf.putAtt(ncid,varid(end),'coordinates','lon lat tair_time');

        %-----------------------------------------------------------------%  
        case 'Qair'
          varid(end+1) = netcdf.defVar(ncid,'qair_time','double',dimid(i+2));
            netcdf.putAtt(ncid,varid(end),'long_name','surface air relative humidity time');
            netcdf.putAtt(ncid,varid(end),'units','days since initialization');
          varid(end+1) = netcdf.defVar(ncid,'Qair','double',[dimid(1) dimid(2) dimid(i+2)]);
            netcdf.putAtt(ncid,varid(end),'long_name','surface air relative humidity');
            netcdf.putAtt(ncid,varid(end),'units','percentage');
            netcdf.putAtt(ncid,varid(end),'time','qair_time');
            netcdf.putAtt(ncid,varid(end),'coordinates','lon lat qair_time'); 

        %-----------------------------------------------------------------%  
        case 'Cfra'
          varid(end+1) = netcdf.defVar(ncid,'cloud_time','double',dimid(i+2));
            netcdf.putAtt(ncid,varid(end),'long_name','cloud fraction time');
            netcdf.putAtt(ncid,varid(end),'units','days since initialization');
          varid(end+1) = netcdf.defVar(ncid,'Cfra','double',[dimid(1) dimid(2) dimid(i+2)]);
            netcdf.putAtt(ncid,varid(end),'long_name','cloud fraction');
            netcdf.putAtt(ncid,varid(end),'units','fraction');
            netcdf.putAtt(ncid,varid(end),'time','cloud_time');
            netcdf.putAtt(ncid,varid(end),'coordinates','lon lat cloud_time');

        %-----------------------------------------------------------------%  
        case 'rain'
          varid(end+1) = netcdf.defVar(ncid,'rain_time','double',dimid(i+2));
            netcdf.putAtt(ncid,varid(end),'long_name','rain fall time');
            netcdf.putAtt(ncid,varid(end),'units','days since initialization');
          varid(end+1) = netcdf.defVar(ncid,'rain','double',[dimid(1) dimid(2) dimid(i+2)]);
            netcdf.putAtt(ncid,varid(end),'long_name','rain fall');
            netcdf.putAtt(ncid,varid(end),'units','kilogram meter-2 second-1');
            netcdf.putAtt(ncid,varid(end),'time','rain_time');
            netcdf.putAtt(ncid,varid(end),'coordinates','lon lat rain_time');  

        %-----------------------------------------------------------------%      
        case {'lwrad','lwrad_down'}
          varid(end+1) = netcdf.defVar(ncid,'lrf_time','double',dimid(i+2));
            netcdf.putAtt(ncid,varid(end),'long_name','longwave radiation flux time');
            netcdf.putAtt(ncid,varid(end),'units','days since initialization');
          varid(end+1) = netcdf.defVar(ncid,var_to_incl{i},'double',[dimid(1) dimid(2) dimid(i+2)]);
            netcdf.putAtt(ncid,varid(end),'long_name','longwave radiation flux');
            netcdf.putAtt(ncid,varid(end),'units','Watts meter-2');
            netcdf.putAtt(ncid,varid(end),'time','lrf_time');
            netcdf.putAtt(ncid,varid(end),'coordinates','lon lat lrf_time');
            netcdf.putAtt(ncid,varid(end),'positive value','downward flux, heating');
            netcdf.putAtt(ncid,varid(end),'negative value','upward flux, cooling');

        %-----------------------------------------------------------------%      
        case {'swrad','swrad_down'}
          varid(end+1) = netcdf.defVar(ncid,'srf_time','double',dimid(i+2));
            netcdf.putAtt(ncid,varid(end),'long_name','shortwave radiation flux time');
            netcdf.putAtt(ncid,varid(end),'units','days since initialization');
          varid(end+1) = netcdf.defVar(ncid,var_to_incl{i},'double',[dimid(1) dimid(2) dimid(i+2)]);
            netcdf.putAtt(ncid,varid(end),'long_name','shortwave radiation flux');
            netcdf.putAtt(ncid,varid(end),'units','Watts meter-2');
            netcdf.putAtt(ncid,varid(end),'time','srf_time');
            netcdf.putAtt(ncid,varid(end),'coordinates','lon lat srf_time');
            netcdf.putAtt(ncid,varid(end),'positive value','downward flux, heating');
            netcdf.putAtt(ncid,varid(end),'negative value','upward flux, cooling');    
        
        %-----------------------------------------------------------------%      
        otherwise
          error(['Unknown or unimplemented variable name: ' var_to_incl{i}]);
    end
    end
    clear i;
   
    % End definitions and close
    netcdf.endDef(ncid);
    netcdf.close(ncid);
  
    % Success!
    status = 1;
  
end