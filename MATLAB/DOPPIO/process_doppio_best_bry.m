%=========================================================================%
% gen_bry_roms_doppio.m
% Generate a ROMS-format boundary forcing file using DOPPIO data
% by Arin Nelson
% on 07/18/2021
% 
% last edited by Arin Nelson on 08/03/2021
%=========================================================================%
clc; clear mex; addpath('../Utilities');

% Time options
date_start = [2019,01,01];	% Start year, month, day to gather data for
date_end   = [2019,12,31];	% End   year, month, day to gather data for

% Grid options
save_gridfile = 'C:/Library/ROMS_Stuff/Resources/ngbay_grd.nc';                                   % ROMS grid file
save_infofile = 'C:/Library/ROMS_Stuff/Resources/forcefiles_riroms/riroms_bry_2006_dt06hrs.nc';   % Contains vertical grid info

% Other options
smooth_factor = 0;      % Low-pass filter with this window (in days).   e.g., to de-tide NECOFS, Dave Ullman uses 36hrs
save_file     = 'F:/OSOM_Data_Repo/_finished/OSOM_bry_2018_DOPPIO_best_old.nc';
save_title    = 'BOUNDARY forcing file: DOPPIO (best) source';
save_WSEN     = [1 1 1 0];      % Include west, south, east, and/or north boundaries

% Variable info (Var Label ,dimensions (2 or 3), grid type 
var_to_get = {'zeta',  2,  1;  ...
              'ubar',  2,  2;  ...
              'vbar',  2,  3;  ...
              'u',     3,  2;  ...
              'v',     3,  3;  ...
              'temp',  3,  1;  ...
              'salt',  3,  1;  ...
             };
n_var = size(var_to_get,1);

% Source info
src_name     = 'DOPPIO_best';
src_dir      = 'F:/OSOM_Data_Repo/DOPPIO/best/';
src_gridfile = 'F:/OSOM_Data_Repo/DOPPIO/doppio_grid.nc';
src_dt       = 1/24;

%=========================================================================%

% OSOM grid variables
src_grid  = grid_get(src_gridfile);
save_grid = grid_get(save_gridfile,save_infofile);

% Grid variables (grid type, boundary)
save_lon   = cell(3,4);     src_lon   = cell(3,4);
save_lat   = cell(3,4);     src_lat   = cell(3,4);
save_mask  = cell(3,4);     src_mask  = cell(3,4);
save_angle = cell(3,4);     src_angle = cell(3,4);

% Grid types and boundaries strings
str_var  = {'lon','lat','mask','angle'};
str_grid = {'rho','u','v'};
str_brdr = {'(1,:)''','(:,1)','(end,:)''','(:,end)'};
for i=1:numel(str_grid)
for j=1:numel(str_brdr)
for k=1:numel(str_var)    
    eval(['save_' str_var{k} '{i,j} = save_grid.' str_var{k} '_' str_grid{i} str_brdr{j} ';']);
end
end
end
clear i j k;

%-------------------------------------------------------------------------%

% Generate boundary file if not yet done so
if(exist(save_file,'file')~=2)
    
    % Generate the file
    opts = struct;
    opts.FileName = save_file;
    opts.WSEN     = save_WSEN;
    opts.Title    = save_title;
    nc_gen_bry_roms(opts,size(save_grid.lon_rho,1),size(save_grid.lon_rho,2),numel(save_grid.s_rho),var_to_get(:,1));

    % Save grid variables
    str_WSEN = {'west','south','east','north'};
    for i=1:2 %numel(str_var)
    for j=1:numel(str_grid)
    for k=1:numel(str_WSEN)   
    if(save_WSEN(k)==1)    
        eval(['ncwrite(save_file,''' str_var{i} '_' str_grid{j} '_' str_WSEN{k} ''',save_' str_var{i} '{j,k});']);
    end
    end
    end
    end
    clear i j k;
    
    % Save vertical variables
    var_to_transfer = {'spherical','Vtransform','Vstretching','theta_s',...
                       'theta_b','Tcline','hc','s_rho','s_w','Cs_r','Cs_w'};
    for i=1:numel(var_to_transfer)
        ncwrite(save_file,var_to_transfer{i},ncread(save_infofile,var_to_transfer{i}));
    end
    clear i;
                   
end

%-------------------------------------------------------------------------%

% Get current timestep
test_time = ncread(save_file,'bry_time');
test_time(test_time > 1e36) = [];
if(isempty(test_time))
    t_on  = datenum(date_start);
    nt_on = 1;
else
    t_on  = datenum(date_start) + test_time(end) + src_dt;
    nt_on = numel(test_time)+1;
end

% Loop through times
while( t_on <= datenum(date_end) )
clc; disp(['On year/month ' num2str(year(t_on)) '/' sprintf('%0.2d',month(t_on)) '...']);

    % Time variables
    year_on  = num2str( year(t_on) );
    month_on = sprintf('%0.2d',month(t_on));
    n_days   = eomday(year(t_on),month(t_on));
    n_times  = n_days / src_dt;
    
    % Loop through variables
    for iv=1:size(var_to_get,1)
    disp([' On variable ' var_to_get{iv,1} '...']); 
    
        % File name for this variable
        src_file = [src_dir '/' var_to_get{iv,1} '/' src_name '_' var_to_get{iv,1} '_' year_on '_' month_on '.nc'];
        
        % Load data
        switch var_to_get{iv,1}
        
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . %
            case 'zeta'
            
                % Test progress
                try        test = ncread(save_file,['zeta_' str_WSEN{find(save_WSEN==1,1,'last')}],[1 nt_on+n_times-1],[1 1]);
                catch err; test=1e37; end
                if(test>1e36)
                
                    % Read in data
                    src_data = ncread(src_file,'zeta');

                    % Loop through times
                    for it=1:size(src_data,3)

                        % Generate interpolant
                        xx = src_grid.lon_rho(src_grid.mask_rho==1);
                        yy = src_grid.lat_rho(src_grid.mask_rho==1);
                        zz = src_data(:,:,it);  zz = zz(src_grid.mask_rho==1);
                        ntrplnt = scatteredInterpolant(xx(:),yy(:),zz(:),'linear');

                        % Loop through borders
                        for ib=1:4
                        if(save_WSEN(ib)==1)
                           save_data = ntrplnt(save_lon{1,ib},save_lat{1,ib});
                           ncwrite(save_file,['zeta_' str_WSEN{ib}],save_data,[1 nt_on+it-1]);
                        end
                        end
                        clear ib;

                        % Clean-up
                        clear xx yy zz ntrplnt;

                    end
                    clear it;

                    % Clean-up
                    clear src_data;
                
                end
                clear test;
                
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . %
            case 'ubar'
                
                % Test progress
                try test = ncread(save_file,['ubar_' str_WSEN{find(save_WSEN==1,1,'last')}],[1 nt_on+n_times-1],[1 1]);
                catch err; test=1e37; end
                if(test>1e36)
                
                    % Read in data
                    src_fil2 = [src_dir '/vbar/' src_name '_vbar_' year_on '_' month_on '.nc'];
                    src_uonu = ncread(src_file,'ubar');
                    src_vonv = ncread(src_fil2,'vbar');

                    % Loop through times
                    imu = find(src_grid.mask_u==1);
                    imv = find(src_grid.mask_v==1);
                    for it=1:size(src_uonu,3)

                        % Generate u interpolant
                        xu  = src_grid.lon_u(imu);
                        yu  = src_grid.lat_u(imu);
                        uonu = src_uonu(:,:,it);    uonu = uonu(imu);
                        ntrplnt_u = scatteredInterpolant(xu(:),yu(:),uonu(:),'linear');

                        % Generate v interpolant
                        xv  = src_grid.lon_v(imv);
                        yv  = src_grid.lat_v(imv);
                        vonv = src_vonv(:,:,it);    vonv = vonv(imv);
                        ntrplnt_v = scatteredInterpolant(xv(:),yv(:),vonv(:),'linear');   

                        % Interp to each other's grids
                        uonv = ntrplnt_u(xv,yv);
                        vonu = ntrplnt_v(xu,yu);

                        % Create z varaibles
                        zonu = uonu + sqrt(-1).*vonu;
                        zonv = uonv + sqrt(-1).*vonv;

                        % Rotate to eastward-northward
                        zonu = zonu .* exp(sqrt(-1).*(-src_grid.angle_u(imu)));
                        zonv = zonv .* exp(sqrt(-1).*(-src_grid.angle_v(imv)));

                        % Generate interpolants
                        ntrplnt_u = scatteredInterpolant(xu(:),yu(:),zonu(:),'linear');
                        ntrplnt_v = scatteredInterpolant(xv(:),yv(:),zonv(:),'linear');

                        % Loop through borders
                        for ib=1:4
                        if(save_WSEN(ib)==1)

                           % Interpolate to border
                           save_u = ntrplnt_u(save_lon{2,ib},save_lat{2,ib});
                           save_v = ntrplnt_v(save_lon{3,ib},save_lat{3,ib});

                           % Rotate
                           save_u = save_u .* exp(sqrt(-1).*(save_angle{2,ib}));
                           save_v = save_v .* exp(sqrt(-1).*(save_angle{3,ib}));

                           % Write to file
                           ncwrite(save_file,['ubar_' str_WSEN{ib}],real(save_u),[1 nt_on+it-1]);
                           ncwrite(save_file,['vbar_' str_WSEN{ib}],imag(save_v),[1 nt_on+it-1]);

                           % Clean-up
                           clear save_u save_v;

                        end
                        end
                        clear ib;

                        % Clean-up
                        clear xu yu uonu ntrplnt_u vonu zonu ...
                              xv yv vonu ntrplnt_v vonv zonv

                    end
                    clear it;

                    % Clean-up
                    clear src_data;
                
                end
                clear test
                
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . %
            case 'u'
                
                % Test progress
                try test = ncread(save_file,['u_' str_WSEN{find(save_WSEN==1,1,'last')}],[1 numel(save_grid.s_rho) nt_on+n_times-1],[1 1 1]);
                catch err; test=1e37; end
                if(test>1e36)
                
                    % Read in data
                    src_fil2 = [src_dir '/v/' src_name '_v_' year_on '_' month_on '.nc'];
                    src_uonu = ncread(src_file,'u');
                    src_vonv = ncread(src_fil2,'v');

                    % Loop through times
                    imu = find(src_grid.mask_u==1);
                    imv = find(src_grid.mask_v==1);
                    for it=1:size(src_uonu,4)

                        % Data at this time
                        uonu = src_uonu(:,:,:,it);
                        vonv = src_vonv(:,:,:,it);

                        % Interpolate in depth
                        uonu_ntrpd = NaN(size(uonu,1),size(uonu,2),numel(save_grid.s_rho));
                        vonv_ntrpd = NaN(size(vonv,1),size(vonv,2),numel(save_grid.s_rho));
                        for ix=1:size(vonv,1)
                        for iy=1:size(uonu,2)

                            % for u
                            if(ix<=size(uonu,1))
                            if(src_grid.mask_u(ix,iy)==1)  
                                uonu_ntrpd(ix,iy,:) = interp1( src_grid.s_rho, squeeze(uonu(ix,iy,:)), save_grid.s_rho, 'linear' ); 
                            end
                            end

                             % for v
                            if(iy<=size(vonv,2))
                            if(src_grid.mask_v(ix,iy)==1)  
                                vonv_ntrpd(ix,iy,:) = interp1( src_grid.s_rho, squeeze(vonv(ix,iy,:)), save_grid.s_rho, 'linear' ); 
                            end
                            end

                        end
                        end
                        clear ix iy uonu vonv;

                        % Save at each depth
                        xu = src_grid.lon_u(imu);
                        yu = src_grid.lat_u(imu);
                        xv = src_grid.lon_v(imv);
                        yv = src_grid.lat_v(imv);
                        for iz=1:numel(save_grid.s_rho)

                            % Data at this depth
                            uonu = vonv_ntrpd(:,:,iz);  uonu = uonu(src_grid.mask_u==1);
                            vonv = vonv_ntrpd(:,:,iz);  vonv = vonv(src_grid.mask_v==1);

                            % Generate interpolants
                            ntrplnt_u = scatteredInterpolant(xu(:),yu(:),uonu(:),'linear');
                            ntrplnt_v = scatteredInterpolant(xv(:),yv(:),vonv(:),'linear');

                            % Interpolate to each others grids
                            uonv = ntrplnt_u(xv,yv);
                            vonu = ntrplnt_v(xu,yu);

                            % Construct complex velocities
                            zonu = uonu + sqrt(-1).*vonu;
                            zonv = uonv + sqrt(-1).*vonv;

                            % Rotate to eastward/northward
                            zonu = zonu .* exp(sqrt(-1).*(-src_grid.angle_u(imu)));
                            zonv = zonv .* exp(sqrt(-1).*(-src_grid.angle_v(imv)));

                            % Generate interpolants
                            ntrplnt_u = scatteredInterpolant(xu(:),yu(:),zonu(:),'linear');
                            ntrplnt_v = scatteredInterpolant(xv(:),yv(:),zonv(:),'linear');

                            % Loop through borders
                            for ib=1:4
                            if(save_WSEN(ib)==1)

                                % Interpolate to border
                                save_u = ntrplnt_u(save_lon{2,ib},save_lat{2,ib});
                                save_v = ntrplnt_v(save_lon{3,ib},save_lat{3,ib});

                                % Rotate
                                save_u = save_u .* exp(sqrt(-1).*(save_angle{2,ib}));
                                save_v = save_v .* exp(sqrt(-1).*(save_angle{3,ib}));

                                % Write to file
                                ncwrite(save_file,['u_' str_WSEN{ib}],real(save_u),[1 1 nt_on+it-1]);
                                ncwrite(save_file,['v_' str_WSEN{ib}],imag(save_v),[1 1 nt_on+it-1]);

                                % Clean-up
                                clear save_u save_v;

                            end
                            end
                            clear ib;

                            % Clean-up
                            clear zz ntrplnt;

                        end
                        clear iz src_ntrpd;

                    end
                    clear it src_data;
                    
                end
                clear test;
                
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . %
            case {'temp','salt'}
                
                % Test progress
                try 
                test = ncread(save_file,[var_to_get{iv,1} '_' str_WSEN{find(save_WSEN==1,1,'last')}],[1 numel(save_grid.s_rho) nt_on+n_times-1],[1 1 1]);
                catch err; test=1e37; end
                if(test>1e36)
                
                    % Read in data
                    src_data = ncread(src_file,var_to_get{iv,1});

                    % Loop through times
                    for it=1:size(src_data,4)

                        % Data at this time
                        z = src_data(:,:,:,it);

                        % Interpolate in depth
                        src_ntrpd = NaN(size(src_data,1),size(src_data,2),numel(save_grid.s_rho));
                        for ix=1:size(src_data,1)
                        for iy=1:size(src_data,2)
                        if(src_grid.mask_rho(ix,iy)==1)
                           src_ntrpd(ix,iy,:) = interp1( src_grid.s_rho, squeeze(z(ix,iy,:)), save_grid.s_rho, 'linear' ); 
                        end
                        end
                        end
                        clear ix iy z;

                        % Save at each depth
                        xx = src_grid.lon_rho(src_grid.mask_rho==1);
                        yy = src_grid.lat_rho(src_grid.mask_rho==1);
                        for iz=1:numel(save_grid.s_rho)

                            % Generate interpolant
                            zz = src_ntrpd(:,:,iz);  zz = zz(src_grid.mask_rho==1);
                            ntrplnt = scatteredInterpolant(xx(:),yy(:),zz(:),'linear');

                            % Loop through borders
                            for ib=1:4
                            if(save_WSEN(ib)==1)
                                save_data = ntrplnt(save_lon{1,ib},save_lat{1,ib});
                                ncwrite(save_file,[var_to_get{iv,1} '_' str_WSEN{ib}],save_data,[1 1 nt_on+it-1]);
                            end
                            end
                            clear ib;

                            % Clean-up
                            clear zz ntrplnt;

                        end
                        clear iz src_ntrpd;

                    end
                    clear it src_data;
                
                end
                
            % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . %   
            otherwise
                % NULL
        end
        clear src_file;
        
    end
    clear iv;
    
    % Once gathered all variables, save time
    t_to_write = (t_on : src_dt : t_on+n_days+1-src_dt) - datenum(date_start);
    ncwrite(save_file,'bry_time',t_to_write,nt_on);
    
    % Onto next month
    t_on = t_on   + n_days + 1;
    nt_on = nt_on + n_times + 1;
    
    % Clean-up
    clear year_on month_on n_days t_to_write;
    
end
clear t_on;

%=========================================================================%