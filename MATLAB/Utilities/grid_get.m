function grid = grid_get(grid_file)

  % Init struct
  grid = struct;

  % Grid variables
  str_var   = {'lon','lat','mask'};
  str_grid  = {'rho','u','v'};
  for i=1:numel(str_var)
  for j=1:numel(str_grid)
    eval(['grid(j).' str_var{i} '=ncread(grid_file,[''' str_var{i} '_' str_grid{j} ''']);']);
  end
  end
  clear i j;
  
  % Angle requires additional calculation
  grid(1).angle = ncread(grid_file,'angle');
  
  % Other 2 angles must be computed
  grid(2).angle = ( grid(1).angle(1:end-1,:) + grid(1).angle(2:end,:) )./2;
  grid(3).angle = ( grid(1).angle(:,1:end-1) + grid(1).angle(:,2:end) )./2;

end