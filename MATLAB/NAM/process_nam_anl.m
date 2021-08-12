%=========================================================================%
% process_nam_anl.m
% Process North American Mesoscale (NAM) model Analysis (ANL) data
% by Arin Nelson
% on 07/31/2021
% 
% Good reference from ROMS variables:
% https://code.usgs.gov/coawstmodel/COAWST/-/blob/84738f369c84b51109c8f9a323c0ddbe9c6f7587/Tools/mfiles/mtools/ncei_2roms.m
% 
% last edited by Arin Nelson on 07/31/2021
%=========================================================================%
clear; clc; addpath('../Utilities');   % clear mex;

% Time Options
date_start = [2018,01,01];	% Start year, month, day to gather data for
date_end   = [2018,12,31];	% End   year, month, day to gather data for
date_plus  = false;         % Have last entry be the first timestep of the day after date_end

% Interpolation options
lon = -72.7 : 0.1 : -69.9;
lat =  40.5 : 0.1 :  42.2;

% Other options
nam_dir       = 'F:/OSOM_Data_Repo/NAM/';
nam_grid      = 'F:/OSOM_Data_Repo/NAM/nam_grid.nc';
var_to_get    = {'wind'};
frcfile_name  = 'OSOM_frc_winds_2018_NAM.nc';
frcfile_title = 'Surface forcing: NAM, spatially variable';
incl_cast     = 1;  % Include forecasted variables


% Shortwave radiation options
swrad_factor    = 1-0.23;             % Multiplicative factor (FROM DOPPIO)
swrad_dailyavg  = false;               % Save swrad as a daily average

% Var info (Var Label , Source , options)
var_info = {'wind',      'wind_time';  ...
            'Tair',      'tair_time';  ...
            'Pair',      'pair_time';  ...
            'Qair',      'qair_time';  ...
            'rain',      'rain_time';  ...
            'lwrad_down','lrf_time';   ...
            'swrad',     'srf_time';   ...
           };
         
%=========================================================================%

% Grid dimensions
nx = numel(lon);
ny = numel(lat);

% Interpolation meshgrid
[latm,lonm] = meshgrid(lat,lon);
    
% Relate variable names to var_info
n_var = numel(var_to_get);
i_var = zeros(n_var,1);
for iv=1:n_var
    i_var(iv) = find( strcmp(var_info(:,1),var_to_get{iv})==1 );
end
clear iv;

% Generate interpolation info it not yet available
if( exist('interp_info_NAM.mat','file')~=2 )

	% Data from grid file  
    nam_lon  = ncread(nam_grid,'lon' );
    nam_lat  = ncread(nam_grid,'lat' );
    nam_mask = ncread(nam_grid,'mask');
    
    % X-Y data
    [    x,    y] = grn2eqa(latm,   lonm   );
    [nam_x,nam_y] = grn2eqa(nam_lat,nam_lon);

    % Inverse-distance-weighted interpolant
    [ntrp_i, ntrp_j, ntrp_w] = grid_interpolant(nam_x,nam_y,x,y);
    
    % Save interp into
    save('interp_info_NAM.mat','x','y','nam_x','nam_y','ntrp_i','ntrp_j','ntrp_w','nam_mask');
          
else
    load('interp_info_NAM.mat');
end

% Input grid size
nam_nx = size(nam_x,1);
nam_ny = size(nam_y,2);

%  Indices of data to get
i0_get = min(ntrp_i(:));     ni_get = max(ntrp_i(:))-i0_get+1;
j0_get = min(ntrp_j(:));     nj_get = max(ntrp_j(:))-j0_get+1;

% Parsed grid
ii = i0_get : (i0_get+ni_get-1);
jj = j0_get : (j0_get+nj_get-1);
nam_x = nam_x(ii,jj);
nam_y = nam_y(ii,jj);
nam_m = nam_mask(ii,jj);
iim   = find(nam_m==1);

% Generate forcing file if it doesn't yet exist
if(exist(frcfile_name,'file')~=2)
    
    % Generate the file
    nc_gen_frc_roms(frcfile_name,nx,ny,var_to_get,frcfile_title);

    % Save grid variables
	ncwrite(frcfile_name,'lon',lon);
    ncwrite(frcfile_name,'lat',lat);
    
end

%-------------------------------------------------------------------------%

    % If interrupted, determine first time step for each variable
    % TO DO

	% Loop through timesteps
	t_on = datenum(date_start);
    it0  = ones(numel(var_to_get),1);
	while(t_on <= datenum(date_end))
    clc; disp(['On date ' datestr(t_on) '...']);
    
        % Year, month, day
        year_on  = num2str( year(t_on) );
        month_on = sprintf( '%0.2d', month(t_on) );
        day_on   = sprintf( '%0.2d', day(t_on) );

        % Loop through wanted variables
        for iv=1:n_var
            
            % Some things are unique to all data files
            nam_file = [nam_dir '/' var_to_get{iv} '/' year_on '/' month_on '/' var_to_get{iv} '_' year_on '_' month_on '_' day_on '.nc'];
            nam_time = ncread(nam_file,'time');
            
            % Determine optimal time order (analysis > forecast)
            avail_times = round(unique(nam_time));     avail_times(avail_times>=24) = [];
            n_times     = numel(avail_times);
            i_times     = cell(n_times,1);
            for i=1:n_times
                [ii,jj] = find( nam_time == avail_times(i) );
                i_times{i} = [ii,jj];
            end
            clear i ii jj;
                    
            % Some other things are variable-dependent
            switch var_to_get{iv}
    
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . %
                case 'wind'

                    % Get data
                    nam_Uwind = ncread(nam_file,'Uwind',[i0_get, j0_get, 1, 1],[ni_get, nj_get, inf, inf]);
                    nam_Vwind = ncread(nam_file,'Vwind',[i0_get, j0_get, 1, 1],[ni_get, nj_get, inf, inf]);
                    
                    % Init outputs
                    out_Uwind = NaN(ni_get,nj_get,0);
                    out_Vwind = NaN(ni_get,nj_get,0);
                    out_Time  = [];
                    
                    % Loop through times
                    for i=1:n_times
                        
                        % Time indices for this time
                        tmp = i_times{i};
                        
                        % Sort in descending order based on columns
                        [~,jj] = sort( tmp(:,1).*100 + tmp(:,2) );
                        tmp = tmp(jj,:);
                        clear jj;
                        
                        % Determine if its valid
                        if(size(tmp,1)~=1)
                            check = false;
                            while(check==false && ~isempty(tmp))
                                testU = nam_Uwind(:,:,tmp(1,1),tmp(1,2));
                                testV = nam_Vwind(:,:,tmp(1,1),tmp(1,2));
                                if( (all(isnan(testU(:))) && all(isnan(testV(:)))) || (all(abs(testU(:))>1e36) && all(abs(testV(:))>1e36)) )
                                    tmp(1,:) = [];
                                else
                                    check = true;
                                end
                                clear testU testV;
                            end
                        end
                        
                        % Data at this time, if available
                        if(~isempty(tmp))
                        	out_Uwind(:,:,end+1) = nam_Uwind(:,:,tmp(1,1),tmp(1,2));
                            out_Vwind(:,:,end+1) = nam_Vwind(:,:,tmp(1,1),tmp(1,2));
                            out_Time(end+1)      = avail_times(i);
                        end
                        clear tmp;
                        
                    end
                    clear i;
                    
                    % Interpolate to output grid
                    ntrp_Uwind = NaN(nx,ny,numel(out_Time));
                    ntrp_Vwind = NaN(nx,ny,numel(out_Time));
                    for ix=1:nx
                    for iy=1:ny
                        wght = squeeze(ntrp_w(ix,iy,:,:));
                        ii   = ntrp_i(ix,iy,:)-i0_get+1;    ii=ii(:);
                        jj   = ntrp_j(ix,iy,:)-j0_get+1;    jj=jj(:);
                        for i=1:numel(out_Time)
                            tmpu = squeeze(out_Uwind(ii,jj,i));
                            tmpv = squeeze(out_Vwind(ii,jj,i));
                            ntrp_Uwind(ix,iy,i) = nansum( tmpu(:).*wght(:) ) ./ nansum( wght(~isnan(tmpu(:))) );
                            ntrp_Vwind(ix,iy,i) = nansum( tmpv(:).*wght(:) ) ./ nansum( wght(~isnan(tmpu(:))) );
                        end
                        clear tmpu tmpv ii jj wgth;
                    end
                    end
                    clear ix iy;
                    
%                     % TEST
%                     ii = i0_get : (i0_get+ni_get-1);
%                     jj = j0_get : (j0_get+nj_get-1);
%                     xx = nam_x(ii,jj);  xLim = [min(xx(:)) max(xx(:))];
%                     yy = nam_y(ii,jj);  yLim = [min(yy(:)) max(yy(:))];
%                     cc = out_Uwind(:,:,1);  cc = [min(cc(:)) max(cc(:))];
%                     subplot(1,2,1); surf(nam_x(ii,jj),nam_y(ii,jj),out_Uwind(:,:,1),'edgecolor','none');  colorbar; view(2);	xlim(xLim); ylim(yLim); caxis(cc);
%                     subplot(1,2,2); surf(x,y,ntrp_Uwind(:,:,1),'edgecolor','none');  colorbar;  view(2); xlim(xLim); ylim(yLim); caxis(cc);
%                     
                    % Save data
                    ncwrite(frcfile_name,'Uwind',    ntrp_Uwind,[1 1 it0(iv)]);
                    ncwrite(frcfile_name,'Vwind',    ntrp_Vwind,[1 1 it0(iv)]);
                    ncwrite(frcfile_name,'wind_time',out_Time,   it0(iv) );
                    it0(iv) = it0(iv) + numel(out_Time);
                    
                    % Clean-up
                    clear nam_Uwind nam_Vwind out_Uwind out_Vwind ntrp_Uwind ntrp_Vwind out_Time;
                 
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . %    
                case {'lwrad','swrad'}
                
                    
                    
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . %    
                case 'rain'
                
                    
                    
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . %    
                case {'Tair','Pair','Qair','cloud'}
                
                    % Get data
                    nam_data = ncread(nam_file,var_to_get{iv},[i0_get, j0_get, 1, 1],[ni_get, nj_get, inf, inf]);

                    % Init outputs
                    out_data = NaN(ni_get,nj_get,0);
                    out_Time = [];
                    
                    % Loop through times
                    for i=1:n_times
                        
                        % Time indices for this time
                        tmp = i_times{i};
                        
                        % Sort in descending order based on columns
                        [~,jj] = sort( tmp(:,1).*100 + tmp(:,2) );
                        tmp = tmp(jj,:);
                        clear jj;
                        
                        % Determine if its valid
                        if(size(tmp,1)~=1)
                            check = false;
                            while(check==false && ~isempty(tmp))
                                test = nam_data(:,:,tmp(1,1),tmp(1,2));
                                if( (all(isnan(testU(:))) && all(isnan(testV(:)))) || (all(abs(testU(:))>1e36) && all(abs(testV(:))>1e36)) )
                                    tmp(1,:) = [];
                                else
                                    check = true;
                                end
                                clear testU testV;
                            end
                        end
                        
                        % Data at this time, if available
                        if(~isempty(tmp))
                        	out_data(:,:,end+1) = nam_data(:,:,tmp(1,1),tmp(1,2));
                            out_Time(end+1)     = avail_times(i);
                        end
                        clear tmp;
                        
                    end
                    clear i;
                    
                    % Interpolate to output grid
                    ntrp_data = NaN(nx,ny,numel(out_Time));
                    for ix=1:nx
                    for iy=1:ny  
                        wght = squeeze(ntrp_w(ix,iy,:,:));
                        ii   = ntrp_i(ix,iy,:)-i0_get+1;    ii=ii(:);
                        jj   = ntrp_j(ix,iy,:)-j0_get+1;    jj=jj(:);
                        wgth = wgth .* nam_m;     % SO INTERPOLATION IS ONLY DONE FROM WATER POINTS!
                        for i=1:numel(out_Time)
                            tmp = squeeze(out_data(ii,jj,i));
                            ntrp_data(ix,iy,i) = nansum( tmpu(:).*wght(:) ) ./ nansum( wght(~isnan(tmpu(:))) );
                        end
                        clear tmpu tmpv ii jj wgth;
                    end
                    end
                    clear ix iy;
                    
                    % Data missed is interpolated from nearest water grid
                    % point
                    for it=1:numel(out_Time)
                       tmp = ntrp_data(:,:,it);
                       dat = out_data(:,:,it);
                       ntrplnt = scatteredInterpolant(nam_x(iim),nam_y(iim),dat(iim),'nearest');
                       tmp(isnan(tmp)) = ntrplnt(x(iim),y(iim));
                    end

                    % Save data
                    ncwrite(frcfile_name,var_to_get{iv},       ntrp_data,[1 1 it0(iv)]);
                    ncwrite(frcfile_name,var_info{i_var(iv),2},out_Time,  it0(iv) );
                    it0(iv) = it0(iv) + numel(out_Time);
                    
                    % Clean-up
                    clear nam_data out_data ntrp_data out_Time;
                    
                    
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . %    
                otherwise
                    error(['Unknown or unimplemented variable name: ' var_to_get{iv}]);
            end
            
            % Clean-up
            clear nam_file nam_time;
            
        end
        clear iv;
        
        % Onto next day
        t_on = t_on + 1;
      
        % Clean-up
        clear year_on month_on day_on;
      
	end
    clear t_on;
            
%-------------------------------------------------------------------------%