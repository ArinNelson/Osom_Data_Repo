function grid = grid_truncate(grid,ii,jj)

  % Variables and grids
  str_var   = {'lon','lat','mask','angle'};
  str_grid  = {'rho','u','v'};
  str_trunc = {'(ii,jj)','(ii(1:end-1),jj)','(ii,jj(1:end-1))'};
  for i=1:numel(str_var)
  for j=1:numel(str_grid)
    eval(['grid(j).' str_var{i} ' = grid(j).' str_var{i} str_trunc{j} ';']);
  end
  end
  clear i j;

end