function grid = grid_border(grid,wsen)

    % truncate string
    str_var   = {'lon','lat','mask','angle'};
    str_grid  = {'rho','u','v'};
    str_trunc = {'(1,:)''','(:,1)','(end,:)''','(:,end)'};

    % Switch based on wsen
    for i=1:numel(str_var)
    for j=1:numel(str_grid)
        eval(['grid.' str_var{i} '_' str_grid{j} '= grid.' str_var{i} '_' str_grid{j} str_trunc{wsen} ';']);
    end
    end

end