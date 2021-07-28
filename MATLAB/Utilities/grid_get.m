function grid = grid_get(grid_file,varargin)

    % Init struct
    grid = struct;

    % if nargin==1, info_file = grid_file
    if(nargin==1);	info_file = grid_file;
    else;           info_file = varargin{1};
    end
  
    % Variable info's
    grid_info = ncinfo(grid_file);    grid_vars = {grid_info.Variables.Name};
    info_info = ncinfo(info_file);    info_vars = {info_info.Variables.Name};
  
    % Horizontal grid variables
    str_var   = {'lon','lat','mask','x','y'};
    str_grid  = {'rho','u','v'};
    for i=1:numel(str_var)
    for j=1:numel(str_grid)
      
        % Read in variable, if it exists
        var_name = [str_var{i} '_' str_grid{j}];
        if( any( strcmp(grid_vars,var_name)==1 ) )  
            eval(['grid.' var_name '=ncread(grid_file,[''' var_name ''']);']);
        end
    
    end
    end
    clear i j;
  
    % Angle requires additional calculation
    try grid.angle_rho = ncread(grid_file,'angle'); catch err; grid.angle_rho = ncread(grid_file,'angle_rho'); end
    grid.angle_u   = ( grid.angle_rho(1:end-1,:) + grid.angle_rho(2:end,:) )./2;
    grid.angle_v   = ( grid.angle_rho(:,1:end-1) + grid.angle_rho(:,2:end) )./2;
    %grid.angle_psi = ( grid.angle_rho(1:end-1,1:end-1) + grid.angle_rho(2:end,2:end) )./2;
    
%     % Other horizontal variables (commented out for now since I haven't had use for them...)
%     var_to_get = {'pm','pn','dndx','dmde','h','f','xl','el','spherical'};
%     for i=1:numel(var_to_get)
%     if( any( strcmp(grid_vars,var_to_get{i})==1 ) )
%         eval(['grid.' var_to_get{i} ' = ncread(grid_file,''' var_to_get{i} ''');']);
%     end
%     end
%     clear i var_to_get;
%   
%     % Compute area
%     grid.area_rho = 1./(grid.pm .* grid.pn);
  
    % Vertical variables are from the 'info' file
    var_to_get = {'s_rho','s_w'}; %,'Vtransform','Vstretching','theta_s','theta_b','Tcline','hc','Cs_r','Cs_w'
    for i=1:numel(var_to_get)
    if( any( strcmp(info_vars,var_to_get{i})==1 ) )
        eval(['grid.' var_to_get{i} ' = ncread(info_file,''' var_to_get{i} ''');']);
    end
    end
    clear i var_to_get;

end