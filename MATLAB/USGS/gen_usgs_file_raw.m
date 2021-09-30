function gen_usgs_file_raw(usgs_id,save_file,varargin)

% % DEBUGGING
% clear mex; close all;
% usgs_id = '01108410';
% save_dir = '.';

%-------------------------------------------------------------------------%
% GATHER SITE METADATA & VARIABLES IN INVENTORY

    % Info of USGS rdb format
    %readtableweb = @(filename)readtable(filename,'Range',end-2:end);
    %webopts = weboptions('ContentReader',readtableweb);
    webopts = weboptions('Timeout',600);

    % Construct URL containing station metadata and inventory
    meta_url = ['https://waterservices.usgs.gov/nwis/site/?format=rdb' ...
                '&sites=' usgs_id ...
                '&siteOutput=expanded&siteStatus=all'];
    meta_rdb = webread(meta_url,webopts);
    meta_rdb = strsplit(meta_rdb,'\n');
    meta_rdb(cellfun(@isempty,meta_rdb)) = [];
    
    % Loop through header lines and keep label -- name pairs
    att_label = {};
    att_name  = {};
    for i=1:numel(meta_rdb)
    if(meta_rdb{i}(1) == '#')
    if(contains(meta_rdb{i},'--'))
        tmp = strsplit(meta_rdb{i},{'--'});
        att_label{end+1} = tmp{1}(2:end);   att_label{end}(att_label{end}==' ')='';
        att_name{end+1}  = tmp{2}(2:end);
    end
    end
    end
    clear i;
    
    % Gather variable values
    n_var     = numel(att_name);
    att_value = cell(n_var,1);
    tmp_label = strsplit(meta_rdb{end-2},'\t','CollapseDelimiters',false);
    tmp_value = strsplit(meta_rdb{end},'\t','CollapseDelimiters',false);
    for i=1:numel(tmp_label)
    if(~isempty(tmp_label{i}))
        att_value{ strcmp(att_label,tmp_label{i})==1 } = tmp_value{i};
    end
    end
    clear i tmp_*;
    
    % If start_date is specified, save it
    if(nargin>2)
        ii = find( strcmp(att_label,'inventory_dt')==1 );
        att_value{ii} = varargin{1};
        clear ii;
    end

    % Remove empty variables
    iempty = find(cellfun(@isempty,att_value)==1);
    att_label(iempty) = [];
    att_name( iempty) = [];
    att_value(iempty) = [];
    clear iempty;

    % Some will be saved as attributes
    att_info = {'agency',                       'agency_cd';        ...
                'site identification number',   'site_no';          ...
                'site name',                    'station_nm';       ...
                'site type',                    'site_tp_cd';       ...
                'site location',                'map_nm';           ...
                'hydrologic unit code',         'huc_cd';           ...
                'established date',             'inventory_dt';     ...
                'timezone',                     'tz_cd';            ...
                'site honors daylight savings?' 'local_time_fg';    ...
                'project number',               'project_no';       ...
               };

    % Others will be saved as variables
    var_info = {'longitude',        'dec_long_va',      'degrees east';         ...
                'latitude',         'dec_lat_va',       'degrees north';        ...
                'altitude',         'alt_va',           'feet above NAVD 1988';	...
                'drainage_area',    'drain_area_va',    'square miles';         ...
               };
           
    % The variable inventory is contained at another url
    inv_url = ['https://waterdata.usgs.gov/nwis/uv?' ...
               'site_no=' usgs_id ...
               '&format=rdb&siteOutput=expanded&siteStatus=all'];
    inv_rdb = webread(inv_url,webopts);
    inv_rdb = strsplit(inv_rdb,'\n');
    
    % The line before the variable inventory begins with '# Data provided for site '
    i_on = 1;
    while( contains(inv_rdb{i_on},'# Data provided for site ')==0 )
        i_on = i_on + 1;
    end
    
    % The next lines give the instrument #, parameter code, and units
    i_on = i_on + 2;
    tmp  = strsplit(inv_rdb{i_on},{'#',' '});
    while( ~all(cellfun(@isempty,tmp)==1) )
        switch tmp{3}
            case '00010';   var_info{end+1,1} = 'temperature';
            case '00060';   var_info{end+1,1} = 'discharge';    
            case '00065';   var_info{end+1,1} = 'gage_height';
            case '00095';   var_info{end+1,1} = 'conductance';
            case '00300';   var_info{end+1,1} = 'dissolved_oxygen';
            case '00400';   var_info{end+1,1} = 'potential_hydrogen';
            otherwise    
                warning(['Unsupported variable label: ' tmp{3} ', skipping...']);
        end
        var_info{end,2} = tmp{2};
        var_info{end,3} = tmp{4};
        for i=5:numel(tmp)
            var_info{end,3} = [var_info{end,3} ' ' tmp{i}];
        end
        i_on = i_on+1;
        tmp  = strsplit(inv_rdb{i_on},{'#',' '});
    end
            
%-------------------------------------------------------------------------%
% CREATE THE NETCDF FILE      
    
    % Init the NetCDF
    if(exist(save_file,'file')==2); delete(save_file);  end
    ncid = netcdf.create(save_file,'NETCDF4');
    
    % NetCDF constants
    nc_global    = netcdf.getConstant('NC_GLOBAL');
    nc_unlimited = netcdf.getConstant('NC_UNLIMITED');
    
    % Define the global attributes
    netcdf.putAtt(ncid,nc_global,'title',  ['USGS streamgage data file (raw), ID # ' usgs_id]);
    netcdf.putAtt(ncid,nc_global,'history',['Created by gen_usgs_file_raw.m on ' datestr(now)]);
    for i=1:size(att_info,1)
        
        % Read in attribute
        ii  = find(strcmp(att_label,att_info{i,2})==1);
        if(~isempty(ii))
            netcdf.putAtt(ncid,nc_global,att_info{i,1},att_value{ii});
        else
            pause(1e-9);
        end
        clear ii;
        
    end
    
    % Define the time dimension
    dimid = netcdf.defDim(ncid,'time',nc_unlimited);
    
    % Single-valued variables
    varid = [];
    for i=1:4
        varid(end+1) = netcdf.defVar(ncid,var_info{i,1},'NC_DOUBLE',[]);
            netcdf.putAtt(ncid,varid(end),'USGS variable',var_info{i,2});
            netcdf.putAtt(ncid,varid(end),'units',var_info{i,3});
            netcdf.defVarDeflate(ncid,varid(end),true,true,2);
            netcdf.defVarFill(ncid,varid(end),false,-999);
    end
    
    % Time variable
    varid(end+1) = netcdf.defVar(ncid,'time','NC_DOUBLE',dimid);
        netcdf.putAtt(ncid,varid(end),'units','matlab datenum');
        netcdf.putAtt(ncid,varid(end),'timezone','UTC');
        netcdf.defVarDeflate(ncid,varid(end),true,true,2);
        netcdf.defVarFill(ncid,varid(end),false,-999);
    
    % Time-dependent variables
    for i=5:size(var_info,1)
        varid(end+1) = netcdf.defVar(ncid,var_info{i,1},'NC_DOUBLE',dimid);
            netcdf.putAtt(ncid,varid(end),'USGS instrument number',var_info{i,2});
            netcdf.putAtt(ncid,varid(end),'units',var_info{i,3});
            netcdf.defVarDeflate(ncid,varid(end),true,true,2);
            netcdf.defVarFill(ncid,varid(end),false,-999);
        varid(end+1) = netcdf.defVar(ncid,[var_info{i,1} '_qcf'],'char',dimid);   
        	netcdf.putAtt(ncid,varid(end),'long_name','quality control flag');
            netcdf.putAtt(ncid,varid(end),'flag_values','e < > R A P');
            netcdf.putAtt(ncid,varid(end),'flag_meanings','see https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-value-qualification-code-uv_rmk_cd');
            netcdf.defVarFill(ncid,varid(end),false,' ');
    end
    
    % End definitions
    netcdf.endDef(ncid);
    netcdf.close(ncid);
    
    % Write single-value variables to new netcdf file
    for i=1:4
        ii = find(strcmp(att_label,var_info{i,2})==1);
        if(~isempty(ii))
            ncwrite(save_file,var_info{i,1},str2double(att_value{ii}));
        end
    end
        
end