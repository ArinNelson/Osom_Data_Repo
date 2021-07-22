function status = nc_gen_nam_data(file_name,ni,nj,var_name)

    % constants
    c_global    = netcdf.getConstant('NC_GLOBAL');
    c_unlimited = netcdf.getConstant('NC_UNLIMITED');

	% NetCDF creation mode options
	ncopts = netcdf.getConstant('NETCDF4');
            
    % Create the NetCDF   
    ncid = netcdf.create(file_name,ncopts);
  
    % Global attributes
    netcdf.putAtt(ncid,c_global,'type','DATA file');
    netcdf.putAtt(ncid,c_global,'source','North American Mesocale (NAM) model Analysis (ANL) product');
    netcdf.putAtt(ncid,c_global,'history',['generated by nc_gen_nam_data.m on ' datestr(now)]);
    netcdf.putAtt(ncid,c_global,'grid','218');

    % Dimensions
    dimid = [];
    dimid(end+1) = netcdf.defDim(ncid,'ni',ni);
    dimid(end+1) = netcdf.defDim(ncid,'nj',nj);
    dimid(end+1) = netcdf.defDim(ncid,'nt',c_unlimited);

    % Time variable
    varid = [];
    varid(end+1) = netcdf.defVar(ncid,'time','double',dimid(3));
        netcdf.putAtt(ncid,varid(end),'units','days since the beginning of this month'); 
    
    % Other variables
    for iv=1:numel(var_name)
    switch var_name{iv}
      
        case 'Uwind'
            varid(end+1) = netcdf.defVar(ncid,'Uwind','double',[dimid(1) dimid(2) dimid(3)]);
            	netcdf.putAtt(ncid,varid(end),'long_name','10m eastward wind velocity');
                netcdf.putAtt(ncid,varid(end),'units','meters second-1');
                
        case 'Vwind'        
            varid(end+1) = netcdf.defVar(ncid,'Vwind','double',[dimid(1) dimid(2) dimid(3)]);
                netcdf.putAtt(ncid,varid(end),'long_name','10m northward wind velocity');
                netcdf.putAtt(ncid,varid(end),'units','meters second-1');
          
        case 'Pair'
            varid(end+1) = netcdf.defVar(ncid,'Pair','double',[dimid(1) dimid(2) dimid(3)]);
                netcdf.putAtt(ncid,varid(end),'long_name','surface air pressure reduced to MSL');
                netcdf.putAtt(ncid,varid(end),'units','mbar');
            
        case 'Tair'
            varid(end+1) = netcdf.defVar(ncid,'Tair','double',[dimid(1) dimid(2) dimid(3)]);
                netcdf.putAtt(ncid,varid(end),'long_name','10m air temperature');
                netcdf.putAtt(ncid,varid(end),'units','degrees Celsius');        
          
        case 'Qair'
            varid(end+1) = netcdf.defVar(ncid,'Qair','double',[dimid(1) dimid(2) dimid(3)]);
                netcdf.putAtt(ncid,varid(end),'long_name','10m air relative humidity');
                netcdf.putAtt(ncid,varid(end),'units','percentage');  
          
        case 'Cfra'
            varid(end+1) = netcdf.defVar(ncid,'Cfra','double',[dimid(1) dimid(2) dimid(3)]);
                netcdf.putAtt(ncid,varid(end),'long_name','cloud fraction');
                netcdf.putAtt(ncid,varid(end),'units','fraction');  
          
        case 'rain'
            varid(end+1) = netcdf.defVar(ncid,'rain','double',[dimid(1) dimid(2) dimid(3)]);
                netcdf.putAtt(ncid,varid(end),'long_name','rainfall rate');
                netcdf.putAtt(ncid,varid(end),'units','kilogram meters+2 seconds-1');  
          
        case 'lwrad_down'
            varid(end+1) = netcdf.defVar(ncid,'lwrad_down','double',[dimid(1) dimid(2) dimid(3)]);
                netcdf.putAtt(ncid,varid(end),'long_name','surface downward longwave radiation');
                netcdf.putAtt(ncid,varid(end),'units','Watts meters+2');
                
        case 'lwrad_up'
            varid(end+1) = netcdf.defVar(ncid,'lwrad_up','double',[dimid(1) dimid(2) dimid(3)]);
                netcdf.putAtt(ncid,varid(end),'long_name','surface upward longwave radiation');
                netcdf.putAtt(ncid,varid(end),'units','Watts meters+2');
          
        case 'swrad_down'
            varid(end+1) = netcdf.defVar(ncid,'swrad_down','double',[dimid(1) dimid(2) dimid(3)]);
                netcdf.putAtt(ncid,varid(end),'long_name','surface downward shortwave radiation');
                netcdf.putAtt(ncid,varid(end),'units','Watts meters+2');
                
        case 'swrad_up'
            varid(end+1) = netcdf.defVar(ncid,'swrad_up','double',[dimid(1) dimid(2) dimid(3)]);
                netcdf.putAtt(ncid,varid(end),'long_name','surface upward shortwave radiation');
                netcdf.putAtt(ncid,varid(end),'units','Watts meters+2');
          
        otherwise
            error(['Unknown or unimplemented variable: ' var_name{iv}]);
            
    end  
    end
    clear iv;
 
    % End definitions and close
    netcdf.endDef(ncid);
    netcdf.close(ncid);
  
	% Success!
    status = 1;

end