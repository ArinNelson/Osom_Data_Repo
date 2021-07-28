function grid = grid_truncate(grid,ii,jj,kk)

    % Arrays
    if(isempty(ii));    ii = 1:size(grid.lon_rho,1);    end
    if(isempty(jj));    jj = 1:size(grid.lon_rho,2);    end
    if(isempty(kk));    kk = 1:numel(grid.s_w);       end
    
    % Grid variables
    var_name = {'lon','lat','mask','x','y','angle','area'};
    grd_name = {'rho','u','v','psi'};
    grd_ii   = {'(ii,jj)','(ii(1:end-1),jj)','(ii,jj(1:end-1))','(ii(1:end-1),jj(1:end-1))'};
    for i=1:numel(var_name)
    for j=1:numel(grd_name)
        
        % this var
        this_var = [var_name{i} '_' grd_name{j}];
        
        % if it exists, truncate it
        if( isfield(grid,this_var)==1 )
            eval(['grid.' this_var '=grid.' this_var grd_ii{j} ';']);
        end
        
    end
    end
    clear i j;
    
    % Other variables (will implement more as they're needed
    var_info = {'pm',       '(ii,jj)';          ...
                'pn',       '(ii,jj)';          ...
                'h',        '(ii,jj)';          ...
                's_rho',    '(kk(1:end-1))';    ...
                's_w',      '(kk)';             ...
                'Cs_r',     '(kk(1:end-1))';    ...
                'Cs_w',     '(kk)';             ...
               };
           
    % If they exist, truncate them       
    for i=1:size(var_info,1)
    if( isfield(grid,var_info{i,1})==1 )
    	eval(['grid.' var_info{i,1} '=grid.' var_info{i,1} var_info{i,2} ';']);
    end
    end
    clear i j;

end