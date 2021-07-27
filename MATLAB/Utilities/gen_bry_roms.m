%=========================================================================%
% gen_bry_roms.m
% Generate a ROMS-format boundary forcing file using the info specified.
% on 07/18/2021
% 
% last edited by Arin Nelson on 07/18/2021
%=========================================================================%
clc; clear mex;

% Switches
Switch    = zeros(9,1);
Switch(1) = 0;      % Compute interpolation info things
Switch(2) = 1;      % Init the boundary file
Switch(3) = 0;      % Gather & interpolate data to boundaries

% Time Options
date_start = [2018,01,01];	% Start year, month, day to gather data for
date_end   = [2020,12,31];	% End   year, month, day to gather data for
date_plus  = false;         % Have last entry be the first timestep of the day after date_end

% Grid options
output_gridfile = 'C:/Library/ROMS_Stuff/Resources/ngbay_grd.nc';                                   % ROMS grid file
output_infofile = 'C:\Library\ROMS_Stuff\Resources\forcefiles_riroms\riroms_bry_2006_dt06hrs.nc';   % Contains vertical grid info

% Other options
smooth_factor = 0;      % Low-pass filter with this window (in days).   e.g., to de-tide NECOFS, Dave Ullman uses 36hrs
frc_fileopts = struct;                         % Options for save file
frc_fileopts.FileName = 'OSOM_bry_2018-20_DOPPIO.nc';
frc_fileopts.Title    = 'BOUNDARY forcing file: DOPPIO (best) source';
frc_fileopts.WSEN     = [1 1 1 0];      % Include west, south, east, and/or north boundaries

% Var info (Var Label , Source , dimensions (2 or 3))
var_to_get = {'zeta',  'DOPPIO',   2;   ...
              'ubar',  'DOPPIO',   2;   ...
              'vbar',  'DOPPIO',   2;   ...
              'u',     'DOPPIO',   3;   ...
              'v',     'DOPPIO',   3;   ...
              'temp',  'DOPPIO',   3;   ...
              'salt',  'DOPPIO',   3;   ...
             };
         
src_info = {'DOPPIO','../../DOPPIO/best/',  'doppio_grid.nc'; ...
           };
         
%=========================================================================%
if(Switch(01))
if(exist('bry_ntrp_info.mat','file')~=2)
    
    % Unique input grids
    input_src = unique(var_to_get(:,2));
    n_src     = numel(input_src);
    
    % Relate to src_info
    input_ndx  = zeros(n_src,1);
    for is=1:n_src
        input_ndx(is) = find( strcmp( src_info(:,1), input_srcs{is} )==1 );
    end
    
	% Load grids
    output_grid = grid_get(output_gridfile);
    %input_grid  = struct; 
    if( numel(input_srcs)==1 )
        %input_grid = grid_get([ src_info{input_ndx,2} src_info{input_ndx,3} ]); 
        load('F:\OSOM_Data_Repo\git\MATLAB\DOPPIO\dopp_gatherinfo.mat','grid_dopp');
        input_grid = grid_dopp; clear grid_dopp;
    else
        % NOT YET IMPLEMENTED
    end
    
    % For interpolation info calculation, set interior points of mask to 0 (only want info at boundaries)
    for i=1:3;  output_grid(i).mask(2:end-1,2:end-1) = 0;   end
    
    % Grid interpolation info
    %ntrp_info = cell(n_src,3);
    if( numel(input_srcs)==1 )
        
        % Do 
        [ntrp_i, ntrp_j, ntrp_w] = grid_intersect(output_grid,input_grid);  
        
        % Save
        ntrp_info = {ntrp_i, ntrp_j, ntrp_w};
        
    else
        % NOT YET IMPLEMENTED
    end
    
    % Backup info?
    save('bry_ntrp_info.mat','ntrp_info','output_grid');

else;   load('bry_ntrp_info.mat');  
end
end
%=========================================================================%
if(Switch(02))        
%if(exist(frc_fileopts.FileName,'file')~=2)
    
	% Init the boundary file
	nc_gen_bry_roms(frc_fileopts,...
                    size(output_grid(1).lon,1),...
                    size(output_grid(1).lon,2),...
                    15,...
                    var_to_get(:,1));

	% Transfer over vertical grid variables
    var_vert = {'spherical','Vtransform','Vstretching','theta_s','theta_b','Tcline','hc','s_rho','s_w','Cs_r','Cs_w'};
    for i=1:numel(var_vert);    ncwrite(frc_fileopts.FileName,var_vert{i},ncread(output_infofile,var_vert{i}));     end
    
    % Transfer over horizontal grid variables
    str_wsen = {'west','south','east','north'};
    var_horz = {};
    for i=1:4
    if(frc_fileopts.WSEN(i)==1)
        var_horz(end+1:end+6) = {['lon_rho_' str_wsen{i}],...
                                 ['lat_rho_' str_wsen{i}],...
                                 ['lon_u_'   str_wsen{i}],...
                                 ['lat_u_'   str_wsen{i}],...
                                 ['lon_v_'   str_wsen{i}],...
                                 ['lat_v_'   str_wsen{i}],...
                                };
    end
    end
    for i=1:numel(var_horz);    ncwrite(frc_fileopts.FileName,var_horz{i},ncread(output_infofile,var_horz{i}));     end
    
%end   
end
%=========================================================================%
if(Switch(03))  
    
	% Determine start date
    test = ncread(frc_fileopts.FileName,'bry_time');
    if(isempty(test))
        year_on  = date_start(1);
        month_on = date_start(2);
    else
        test     = test(end) + datenum(date_start(1),date_start(2),date_start(3));
        year_on  = year(tmp);
        month_on = month(tmp);
    end
    clear test;
    
    % Loop
    while( datenum(year_on,month_on,1) <= datenum(year_on,month_on,eomday(year_on,month_on)) )
        
        % Loop through wanted variables
        for iv=1:size(var_to_get,1)
            
            % Do transfer in a separate function
            bry_ntrp_fcn(src_info{input_ndx(iv),2}, ...
                         var_to_get{iv},            ...
                         [year_on, month_on],       ...
                         ntrp_info,                 ...
                         frc_fileopts.WSEN,         ...
                         frc_fileopts.FileName      );
            
        end
        clear iv;
        
        % Onto next year/month
        month_on = month_on + 1;
        if(month_on == 13)
            year_on  = year_on + 1;
            month_on = 1;
        end
        
    end
    clear year_on month_on;
    

end
%=========================================================================%      