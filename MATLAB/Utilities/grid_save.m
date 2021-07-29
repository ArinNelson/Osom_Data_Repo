function grid_save(file_name,grid)

    % Dimension sizes
    ni = size(grid.lon_rho,1);
    nj = size(grid.lon_rho,2);
    ns = numel(grid.s_rho);

    % Unlimited dimension
    c_unlimited = netcdf.getConstant('NC_UNLIMITED');
    c_global    = netcdf.getConstant('GLOBAL');

    % Init
    if(exist(file_name,'file')==2); delete(file_name); end
    ncid = netcdf.create(file_name,'NETCDF4');
  
    % Set the global attributes
    netcdf.putAtt(ncid,c_global,'type','GRID file');
    netcdf.putAtt(ncid,c_global,'history',['generated by grid_save.m on ' datestr(now)]);

    % Set the dimensions
    dimid = [];
    dimid(end+1) = netcdf.defDim(ncid,'xi_rho',     ni          );
    dimid(end+1) = netcdf.defDim(ncid,'xi_u',       ni-1        );
    dimid(end+1) = netcdf.defDim(ncid,'xi_v',       ni          );
    dimid(end+1) = netcdf.defDim(ncid,'eta_rho',    nj          );
    dimid(end+1) = netcdf.defDim(ncid,'eta_u',      nj          );
    dimid(end+1) = netcdf.defDim(ncid,'eta_v',      nj-1        );
    dimid(end+1) = netcdf.defDim(ncid,'s_rho',      ns          );
    dimid(end+1) = netcdf.defDim(ncid,'s_w',        ns+1        );

	% Define the vertical grid variables
    varid = [];
    varid(end+1) = netcdf.defVar(ncid,'s_rho','NC_DOUBLE',dimid(end-1));
        netcdf.putAtt(ncid,varid(end),'long_name','S-coordinate at RHO-points');
        netcdf.putAtt(ncid,varid(end),'valid_min','-1');
        netcdf.putAtt(ncid,varid(end),'valid_max','0');
        netcdf.putAtt(ncid,varid(end),'positive','up');
        netcdf.putAtt(ncid,varid(end),'standard_name','ocean_s_coordinate_g2');
        netcdf.putAtt(ncid,varid(end),'formula_terms','s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc');
    varid(end+1) = netcdf.defVar(ncid,'s_w','NC_DOUBLE',dimid(end));
        netcdf.putAtt(ncid,varid(end),'long_name','S-coordinate at W-points');
        netcdf.putAtt(ncid,varid(end),'valid_min','-1');
        netcdf.putAtt(ncid,varid(end),'valid_max','0');
        netcdf.putAtt(ncid,varid(end),'positive','up');
        netcdf.putAtt(ncid,varid(end),'standard_name','ocean_s_coordinate_g2');
        netcdf.putAtt(ncid,varid(end),'formula_terms','s: s_w C: Cs_w eta: zeta depth: h depth_c: hc');

    % Define the horizontal grid variables
    var_prfx = {'lon','lat','mask','angle','x','y'};
    var_grid = {'rho','u','v'};
    var_dims = {dimid([1 4]), dimid([2 5]), dimid([3 6])};
    for i=1:numel(var_prfx)
    for j=1:numel(var_grid)
        
       % Variable
       var_name = [var_prfx{i} '_' var_grid{j}];
       
       % If it exists, define it
       if( isfield(grid,var_name)==1 )
            varid(end+1) = netcdf.defVar(ncid,var_name,'NC_DOUBLE',var_dims{j});
       end
        
        
    end
    end
    clear i j;
    
    % End definitions and close
    netcdf.endDef(ncid);
    netcdf.close(ncid);
    
    % Now write out the variables
    info=ncinfo(file_name);
    for i=1:numel(info.Variables)
        eval(['ncwrite(file_name,''' info.Variables(i).Name ''',grid.' info.Variables(i).Name ');']);
    end
    
    % Same for 
    
end