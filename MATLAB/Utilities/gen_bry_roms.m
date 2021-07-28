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
Switch(2) = 0;      % Init the boundary file
Switch(3) = 1;      % Gather & interpolate data to boundaries

% Time Options
date_start = [2018,01,01];	% Start year, month, day to gather data for
date_end   = [2020,12,31];	% End   year, month, day to gather data for
date_plus  = false;         % Have last entry be the first timestep of the day after date_end

% Grid options
output_gridfile = 'C:/Library/ROMS_Stuff/Resources/ngbay_grd.nc';                                   % ROMS grid file
output_infofile = 'C:/Library/ROMS_Stuff/Resources/forcefiles_riroms/riroms_bry_2006_dt06hrs.nc';   % Contains vertical grid info

% Other options
smooth_factor = 0;      % Low-pass filter with this window (in days).   e.g., to de-tide NECOFS, Dave Ullman uses 36hrs
bry_fileopts = struct;	% Options for save file
bry_fileopts.FileName = 'OSOM_bry_2018-20_DOPPIO.nc';
bry_fileopts.Title    = 'BOUNDARY forcing file: DOPPIO (best) source';
bry_fileopts.WSEN     = [1 1 1 0];      % Include west, south, east, and/or north boundaries

% Variable info (Var Label , Source , dimensions (2 or 3))
var_to_get = {'zeta',  'DOPPIO',   2;   ...
              'ubar',  'DOPPIO',   2;   ...
              'vbar',  'DOPPIO',   2;   ...
              'u',     'DOPPIO',   3;   ...
              'v',     'DOPPIO',   3;   ...
              'temp',  'DOPPIO',   3;   ...
              'salt',  'DOPPIO',   3;   ...
             };
n_var = size(var_to_get,1);

% Source info
src_info = {'DOPPIO_best','F:/OSOM_Data_Repo/DOPPIO/best/',  'F:/OSOM_Data_Repo/DOPPIO/doppio_grid.nc'; ...
           };
         
%=========================================================================%
if(Switch(01))
    
    % Unique input grids
    input_src = unique(var_to_get(:,2));
    n_src     = numel(input_src);
    
    % Relate to src_info
    input_ndx  = zeros(n_src,1);
    for is=1:n_src
        input_ndx(is) = find( strcmp( src_info(:,1), input_src{is} )==1 );
    end
    clear is;
    
	% Load grids
    output_grid = grid_get(output_gridfile,output_infofile);
    input_grid  = struct; 
    for is=1:n_src
        tmp1 = grid_get(src_info{input_ndx(is),3});
        tmp2 = fields(tmp1);
        for j=1:numel(tmp2)
            eval(['input_grid(is).' tmp2{j} '=tmp1.' tmp2{j} ';']);
        end
        clear tmp1 tmp2 j;
    end
    clear is;
    
    % For interpolation info calculation, set interior points of mask to 0 (only want info at boundaries)
    output_grid.mask_rho(2:end-1,2:end-1) = 0;
    output_grid.mask_u(  2:end-1,2:end-1) = 0;
    output_grid.mask_v(  2:end-1,2:end-1) = 0;
    
    % Grid interpolation info
    ntrp_info = cell(n_src,1);
    for is=1:n_src
    if(exist(['bry_ntrp_info_' input_src{is} '.mat'],'file')~=2)
        
        % Gather interpolation info
        [ntrp_i, ntrp_j, ntrp_w] = grid_intersect(output_grid,input_grid(is));  
        
         % Save to variable
         ntrp_info{is} = {ntrp_i, ntrp_j, ntrp_w};
         
         % Also backup to file
         save(['bry_ntrp_info_' input_src{is} '.mat'],'ntrp_info');

    else
        
        % Load previously-defined interpolation info
        tmp = load(['bry_ntrp_info_' input_src{is} '.mat']);  
        ntrp_info{is} = tmp.ntrp_info;
        clear tmp;
        
    end
    end
    clear is;
    
end
%=========================================================================%
if(Switch(02))        
if(exist(bry_fileopts.FileName,'file')~=2)
    
	% Init the boundary file
	nc_gen_bry_roms(bry_fileopts,...
                    size(output_grid.lon_rho,1),...
                    size(output_grid.lon_rho,2),...
                    numel(output_grid.s_rho),...
                    var_to_get(:,1));

	% Transfer over vertical grid variables
    var_vert = {'spherical','Vtransform','Vstretching','theta_s','theta_b','Tcline','hc','s_rho','s_w','Cs_r','Cs_w'};
    for i=1:numel(var_vert);    ncwrite(bry_fileopts.FileName,var_vert{i},ncread(output_infofile,var_vert{i}));     end
    
    % Transfer over horizontal grid variables
    str_wsen = {'west','south','east','north'};
    var_horz = {};
    for i=1:4
    if(bry_fileopts.WSEN(i)==1)
        var_horz(end+1:end+6) = {['lon_rho_' str_wsen{i}],...
                                 ['lat_rho_' str_wsen{i}],...
                                 ['lon_u_'   str_wsen{i}],...
                                 ['lat_u_'   str_wsen{i}],...
                                 ['lon_v_'   str_wsen{i}],...
                                 ['lat_v_'   str_wsen{i}],...
                                };
    end
    end
    for i=1:numel(var_horz);    ncwrite(bry_fileopts.FileName,var_horz{i},ncread(output_infofile,var_horz{i}));     end
    
end
end
%=========================================================================%
if(Switch(03))  
    
	% Determine start date
    test = ncread(bry_fileopts.FileName,'bry_time');
    if(isempty(test))
        year_on  = date_start(1);
        month_on = date_start(2);
    else
        test     = test(end) + datenum(date_start(1),date_start(2),date_start(3));
        year_on  = year(tmp);
        month_on = month(tmp);
    end
    clear test;
    
    % Bounding boxes
    i0 = zeros(4,3,n_src);    ni = zeros(4,3,n_src);
    j0 = zeros(4,3,n_src);    nj = zeros(4,3,n_src);
    for is=1:n_src
    for ig=1:3    
        
        % Interpolation point indices and weights
        ii = ntrp_info{is}{1}{ig};
        jj = ntrp_info{is}{2}{ig};
        %kk = ntrp_info{is}{3}{ig};
        
        % Western edge?
        if(bry_fileopts.WSEN(1)==1)
        	iii = squeeze(ii(:,1,:));   iii = iii(:);   i0(1,ig,is) = min(iii);     ni(1,ig,is) = max(iii)-i0(1,ig,is)+1;
            jjj = squeeze(jj(:,1,:));   jjj = jjj(:);   j0(1,ig,is) = min(jjj);     nj(1,ig,is) = max(jjj)-j0(1,ig,is)+1;
        end
        
        % Southern edge?
        if(bry_fileopts.WSEN(2)==1)
        	iii = squeeze(ii(1,:,:));   iii = iii(:);   i0(2,ig,is) = min(iii);     ni(2,ig,is) = max(iii)-i0(2,ig,is)+1;
            jjj = squeeze(jj(1,:,:));   jjj = jjj(:);   j0(2,ig,is) = min(jjj);     nj(2,ig,is) = max(jjj)-j0(2,ig,is)+1;
        end
        
        % Eastern edge?
        if(bry_fileopts.WSEN(3)==1)
        	iii = squeeze(ii(:,end,:));   iii = iii(:);   i0(3,ig,is) = min(iii);     ni(3,ig,is) = max(iii)-i0(3,ig,is)+1;
            jjj = squeeze(jj(:,end,:));   jjj = jjj(:);   j0(3,ig,is) = min(jjj);     nj(3,ig,is) = max(jjj)-j0(3,ig,is)+1;
        end
        
        % Northern edge?
        if(bry_fileopts.WSEN(4)==1)
        	iii = squeeze(ii(end,:,:));   iii = iii(:);   i0(4,ig,is) = min(iii);     ni(4,ig,is) = max(iii)-i0(4,ig,is)+1;
            jjj = squeeze(jj(end,:,:));   jjj = jjj(:);   j0(4,ig,is) = min(jjj);     nj(4,ig,is) = max(jjj)-j0(4,ig,is)+1;
        end
        
    end
    end
    clear is ig ii jj iii jjj;
    
    % Loop through years & months
    while( datenum(year_on,month_on,1) <= datenum(year_on,month_on,eomday(year_on,month_on)) )
        
        % Loop through variables
        for iv=1:n_var
        
            % This variable's file
            input_filename = [src_info{input_ndx(iv),2} var_info{iv,1} '/' src_info{input_ndx(iv),1} '_' var_info{iv,1} '_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
        
            % If first variable, get time & ensure its within specified time range
            if(iv==1)
                
                % Input time info
                input_time = ncread(input_filename,'time');
                input_it0  = find(time >= datenum(date_start(1),date_start(2),date_start(3)),'first');
                input_itf  = find(time <= datenum(date_end(1),date_end(2),date_end(3))+1,'first');
                input_nt   = itf-it0+1;
                
                % Where to start saving to output file
                output_time = ncread(bry_fileopts.FileName,'bry_time');
                if(isempty(output_time)==1)
                    output_it0 = 1;
                else
                    output_it0 = find(input_time < output_time(1), 1, 'last')+1;
                end
                
            end
            
            %-------------------------------------------------------------%
            switch var_info{iv,1}
                
                case 'zeta'
                    
                    
                case 'ubar'
                    
                
                case 'temp'
                    
                    
                case 'salt'
                    
                    
                case 'u'
           
            
                otherwise
                    %NOTHING
            end
            
            %-------------------------------------------------------------%
            
            % If last variable, save time to bry_time
            if(iv==n_var)
                ncwrite(bry_fileopts.FileName,'bry_time',output_time,output_it0);
            end
        
        end
        
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