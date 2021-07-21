%=========================================================================%
% gen_frc_roms.m
% Generate a ROMS-format surface forcing file using the info specified.
% on 07/18/2021
% 
% last edited by Arin Nelson on 07/18/2021
%=========================================================================%
clc; clear mex;

% Switches
Switch    = zeros(9,1);
Switch(1) = 1;      % Compute interpolation info things
Switch(2) = 1;      % Generate the file & gather the data

% Time Options
date_start = [2010,01,01];	% Start year, month, day to gather data for
date_end   = [2010,12,31];	% End   year, month, day to gather data for
date_plus  = true;          % Have last entry be the first timestep of the day after date_end

% Interpolation options
lon = -74.0 : 0.1 : -69.0;
lat =  40.0 : 0.1 :  42.5;

% Other options
frc_fileopts = struct;                         % Options for save file
frc_fileopts.FileName = 'OSOM_frc_swrad_2010_NAM.nc';                   % Save file name
frc_fileopts.Title    = 'Wind forcing: NAM winds, spatially variable';	% Title attribute of frc file

% do_doppio_fix = true;               % Reduce shortwave radiation by 23%
% nam_grid      = 'nam_grid.nc';      % NAM grid file
% data_dir      = '.\';               % Head directory of NAM data (gathered by gather_nam_anl.m)
% lwrad_type    = 'downward';         % upward, downward, or net longwave radiation
% swrad_type    = 'net';              % upward, downward, or net shortwave radiation
% average_swrad = true;               % Compute daily averages of shortwave radiation
% var_incl      = {'winds'};          % Variables to include in file

% Shortwave radiation options
swrad_factor    = 1-0.23;             % Multiplicative factor (FROM DOPPIO)
swrad_dailyavg  = true;               % Save swrad as a daily average

% Var info
% Var Label , Source , options
% var_to_get = {'winds','NAM',{}; ...
% var_to_get = {'Tair','NAM','tair_time'; ...
%               'Pair','NAM','pair_time'; ...
%               'Qair','NAM','qair_time'; ...
%               'Cfra','NAM','cloud_time'; ...
%               'rain','NAM','rain_time'; ...
%               'lwrad_down','NAM','lrf_time'; ...
var_to_get = {'swrad','NAM','srf_time'; ...
             };
         
%-------------------------------------------------------------------------%

% Constants
% Source name , data directory , grid file name , dt in days
src_info = {'NAM','../../NAM/',  'nam_grid.nc', 0.125;  ...
            'NARR','../../NARR/','narr_grid.nc',0.125; ...
           };

%=========================================================================%
if(Switch(1)==1)

  % Grid dimensions
  nx = numel(lon);
  ny = numel(lat);

  % Interpolation meshgrid
  [latm,lonm] = meshgrid(lat,lon);

  % Load grid info for unique sources
  n_src = size(var_to_get,1);
  src   = struct;
  for i=1:n_src
    
    % This source's entry in src_info
    ii = find( strcmp(src_info(:,1),var_to_get{i,2})==1 );
    if(isempty(ii)); error(['Unknown data source: ' var_to_get{i,2}]); end  
    src(i).Source = src_info{ii,1};
    src(i).Dir    = src_info{ii,2};
    src(i).dt     = src_info{ii,4};
    
    % Data from grid file  
    src(i).x    = ncread([src_info{ii,2} src_info{ii,3}],'x'   );
    src(i).y    = ncread([src_info{ii,2} src_info{ii,3}],'y'   );
    src(i).lon  = ncread([src_info{ii,2} src_info{ii,3}],'lon' );
    src(i).lat  = ncread([src_info{ii,2} src_info{ii,3}],'lat' );
    src(i).m    = ncread([src_info{ii,2} src_info{ii,3}],'mask');
  
    % Source lon/lat meshgrid (if lon & lat are 1D)
    if( any(size(src(i).lon)==1) )
      [tmpY, tmpX] = meshgrid(src(i).lat,src(i).lon);
      src(i).lon   = tmpX;    clear tmpX;
      src(i).lat   = tmpY;    clear tmpY;
    end
  
    % Source x/y 2D grid (if x & y are 1D)
    if( any(size(src(i).x)==1) )
      [tmpY, tmpX] = meshgrid(src(i).y,src(i).x);
      src(i).x     = tmpX;    clear tmpX;
      src(i).y     = tmpY;    clear tmpY;
    end
  
    % Estimate x and y at interpolation points
    ntrplnt_x = scatteredInterpolant(src(i).lon(:),cosd(src(i).lat(:)),src(i).x(:),'linear');
    ntrplnt_y = scatteredInterpolant(src(i).lon(:),cosd(src(i).lat(:)),src(i).y(:),'linear');
    src(i).xi = ntrplnt_x(lonm,cosd(latm));  clear ntrplnt_x;
    src(i).yi = ntrplnt_y(lonm,cosd(latm));  clear ntrplnt_y;
  
  end
  clear i;

  %---INEFFICIENT, BUT NOT TOO SLOW, & BASED ON OLDER CODE---%
  
  % Time array for this file
  tI = datenum(date_start(1),date_start(2),date_start(3),0,0,0);
  tF = datenum(date_end(1),  date_end(2),  date_end(3),  0,0,0);
  if(date_plus==true); tF = tF+1; end
  t  = tI : 1 : tF;  
  
  % Time info
  time_todo = [year(t); month(t);];
  time_todo = unique(time_todo','rows');
  clear tI tF t;

end
%=========================================================================%
if(Switch(2)==1)
  
  % Generate forcing file if it doesn't yet exist
  if(exist(frc_fileopts.FileName,'file')~=2)  
    nc_gen_frc_roms(frc_fileopts,nx,ny,var_to_get(:,1));
    ncwrite(frc_fileopts.FileName,'lon',lon);
    ncwrite(frc_fileopts.FileName,'lat',lat);
  end
  
  % Loop through years and months
  for it=1:size(time_todo,1)
  clc; disp(['On year/month ' num2str(time_todo(it,1)) '/' sprintf('%0.2d',time_todo(it,2)) '...']);
    
    % Loop through wanted variables
    for iv=1:size(var_to_get,1)
    switch var_to_get{iv,1}
    
      %-------------------------------------------------------------------%
      case 'winds'
       
        % Load winds data for this year & month
        input_filename = [src(iv).Dir 'winds/' src(iv).Source '_winds_' num2str(time_todo(it,1)) '_' sprintf('%0.2d',time_todo(it,2)) '.nc'];
        input_Uair = ncread(input_filename,'Uair');
        input_Vair = ncread(input_filename,'Vair');
        input_time = ncread(input_filename,'wind_time') - datenum(date_start(1),date_start(2),date_start(3)) + datenum(time_todo(it,1),time_todo(it,2),1);
        input_dt   = diff(input_time(1:2));
        nii        = numel(input_time);   % size(input_Uair,3)
          
        % Interpolate input to new grid
        Uair = zeros(nx,ny,nii);
        Vair = zeros(nx,ny,nii);
        for i=1:nii
          zz = input_Uair(:,:,i);  
          ntrplnt_U = scatteredInterpolant(src(iv).x(:),src(iv).y(:),zz(:),'linear');
          Uair(:,:,i) = ntrplnt_U(src(iv).xi,src(iv).yi);
          zz = input_Vair(:,:,i);
          ntrplnt_V = scatteredInterpolant(src(iv).x(:),src(iv).y(:),zz(:),'linear');
          Vair(:,:,i) = ntrplnt_V(src(iv).xi,src(iv).yi);
        end
              
        % Before saving, get current time index
        tmp = ncread(frc_fileopts.FileName,'wind_time');
        ii0 = numel(tmp)+1;
          
        % Save to file
        ncwrite(frc_fileopts.FileName,'wind_time',input_time,ii0);
        ncwrite(frc_fileopts.FileName,'Uwind',Uair,[1 1 ii0]);
        ncwrite(frc_fileopts.FileName,'Vwind',Vair,[1 1 ii0]);
          
        % Clean-up
        clear input_* nii ii0 tmp i zz Uair Vair;

      %-------------------------------------------------------------------%  
      case {'Tair','Pair','Qair','Cfra','rain','lwrad_down','swrad_down','lwrad','swrad'}
            
        % Load data for this year & month
        if( strcmp(var_to_get{iv,1},'swrad')==1 || strcmp(var_to_get{iv,1},'lwrad')==1 )
          input_filename = [src(iv).Dir var_to_get{iv,1} '/' src(iv).Source '_' var_to_get{iv,1} '_' num2str(time_todo(it,1)) '_' sprintf('%0.2d',time_todo(it,2)) '.nc'];
          input_data = ncread(input_filename,[var_to_get{iv,1} '_down']) - ncread(input_filename,[var_to_get{iv,1} '_up']);   
        elseif( strcmp(var_to_get{iv,1},'swrad_down')==1 || strcmp(var_to_get{iv,1},'lwrad_down')==1 )
          input_filename = [src(iv).Dir var_to_get{iv,1}(1:5) '/' src(iv).Source '_' var_to_get{iv,1}(1:5) '_' num2str(time_todo(it,1)) '_' sprintf('%0.2d',time_todo(it,2)) '.nc'];
          input_data = ncread(input_filename,var_to_get{iv,1});   
        else
          input_filename = [src(iv).Dir var_to_get{iv,1} '/' src(iv).Source '_' var_to_get{iv,1} '_' num2str(time_todo(it,1)) '_' sprintf('%0.2d',time_todo(it,2)) '.nc'];
          input_data = ncread(input_filename,var_to_get{iv,1}); 
        end
        input_time = ncread(input_filename,'wind_time') - datenum(date_start(1),date_start(2),date_start(3)) + datenum(time_todo(it,1),time_todo(it,2),1);
        input_dt   = diff(input_time(1:2));
        nii        = numel(input_time);
          
        % Interpolate input to new grid (but only from water points!)
        data = zeros(nx,ny,nii);
        imsk = find(src(iv).m == 1);
        for i=1:nii
          zz = input_data(:,:,i);  
          ntrplnt = scatteredInterpolant(src(iv).x(imsk),src(iv).y(imsk),zz(imsk),'linear','nearest');
          data(:,:,i) = ntrplnt(src(iv).xi,src(iv).yi);
        end
        
        % Before saving, get current time index
        tmp = ncread(frc_fileopts.FileName,var_to_get{iv,3});
        ii0 = numel(tmp)+1;
        
        % If swrad, check additional options
        if( strcmp(var_to_get{iv,1},'swrad_down')==1 || strcmp(var_to_get{iv,1},'swrad')==1 )
          data = data.*swrad_factor;
          if swrad_dailyavg==true
            n_avg    = 1/diff(input_time([1 2]));
            data_avg = zeros(nx,ny,nii/n_avg);
            for i=1:(nii/n_avg)
              i_avg = (i-1)*n_avg+1 : i*n_avg;
              data_avg(:,:,i) = nanmean(data(:,:,i_avg),3);
            end
            data = data_avg;
            input_time = input_time(1:n_avg:end);
          end
        end
        
        % Save to file
        ncwrite(frc_fileopts.FileName,var_to_get{iv,3},input_time,ii0);
        ncwrite(frc_fileopts.FileName,var_to_get{iv,1},data,[1 1 ii0]);
        
        % Clean-up
        clear input_* nii data i zz ntrplnt tmp ii0;
          
      %-------------------------------------------------------------------%
      otherwise
          error(['Unknown or unimplemented variable name: ' var_to_get{iv,1}]);  
          
    end
    end
    clear iv;    
    
  end
  clear it;

end
%=========================================================================%