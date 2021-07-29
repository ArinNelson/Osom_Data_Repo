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
var_to_get = {'zeta',  'DOPPIO_best',   2;   ...
              'ubar',  'DOPPIO_best',   2;   ...
              'vbar',  'DOPPIO_best',   2;   ...
              'u',     'DOPPIO_best',   3;   ...
              'v',     'DOPPIO_best',   3;   ...
              'temp',  'DOPPIO_best',   3;   ...
              'salt',  'DOPPIO_best',   3;   ...
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
    test = ncread(bry_fileopts.FileName,'bry_time');    test(test>1e36) = [];
    if(isempty(test))
        year_on  = date_start(1);
        month_on = date_start(2);
    else
        test     = test(end) + datenum(date_start(1),date_start(2),date_start(3));
        year_on  = year(test);
        month_on = month(test);
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
        	iii = squeeze(ii(1,:,:));   iii = iii(:);   i0(1,ig,is) = min(iii);     ni(1,ig,is) = max(iii)-i0(1,ig,is)+1;
            jjj = squeeze(jj(1,:,:));   jjj = jjj(:);   j0(1,ig,is) = min(jjj);     nj(1,ig,is) = max(jjj)-j0(1,ig,is)+1;
        end
        
        % Southern edge?
        if(bry_fileopts.WSEN(2)==1)
        	iii = squeeze(ii(:,1,:));   iii = iii(:);   i0(2,ig,is) = min(iii);     ni(2,ig,is) = max(iii)-i0(2,ig,is)+1;
            jjj = squeeze(jj(:,1,:));   jjj = jjj(:);   j0(2,ig,is) = min(jjj);     nj(2,ig,is) = max(jjj)-j0(2,ig,is)+1;
        end
        
        % Eastern edge?
        if(bry_fileopts.WSEN(3)==1)
        	iii = squeeze(ii(end,:,:));   iii = iii(:);   i0(3,ig,is) = min(iii);     ni(3,ig,is) = max(iii)-i0(3,ig,is)+1;
            jjj = squeeze(jj(end,:,:));   jjj = jjj(:);   j0(3,ig,is) = min(jjj);     nj(3,ig,is) = max(jjj)-j0(3,ig,is)+1;
        end
        
        % Northern edge?
        if(bry_fileopts.WSEN(4)==1)
        	iii = squeeze(ii(:,end,:));   iii = iii(:);   i0(4,ig,is) = min(iii);     ni(4,ig,is) = max(iii)-i0(4,ig,is)+1;
            jjj = squeeze(jj(:,end,:));   jjj = jjj(:);   j0(4,ig,is) = min(jjj);     nj(4,ig,is) = max(jjj)-j0(4,ig,is)+1;
        end
        
    end
    end
    clear is ig ii jj iii jjj;
    
    % Loop through years & months
    while( datenum(year_on,month_on,1) <= datenum(year_on,month_on,eomday(year_on,month_on)) )
        
        % Loop through variables
        for iv=1:n_var
            
            % Variable's src
            i_src = find( strcmp( src_info(:,1), var_to_get{iv,2} )==1 );
        
            % This variable's file
            input_filename = [src_info{i_src,2} var_to_get{iv,1} '/' src_info{i_src,1} '_' var_to_get{iv,1} '_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
        
            % If first variable, get time & ensure its within specified time range
            if(iv==1)
                
                % Input time info
                input_time = ncread(input_filename,'time');
                input_it0  = find(input_time >= datenum(date_start(1),date_start(2),date_start(3)),1,'first');
                input_itf  = find(input_time <= datenum(date_end(1),  date_end(2),  date_end(3))+1,1,'last');
                input_nt   = input_itf-input_it0+1;
                
                % Where to start saving to output file
                output_time = ncread(bry_fileopts.FileName,'bry_time');
                if(isempty(output_time)==1)
                    output_it0 = 1;
                else
                    output_it0 = find(input_time < output_time(1), 1, 'last')+1;
                end
                
            end
            
            switch var_to_get{iv,1}
                
                %% ----------------------------------------------------- %%
                case 'zeta'
                    
                    for i=1:4
                    if(bry_fileopts.WSEN(i)==1)
                        
                        % Gather data surrounding boundary
                        switch i
                            
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            case 1  % West
                                
                                % Input data
                                zeta_in  = ncread(input_filename, 'zeta', [i0(i,1,i_src),j0(i,1,i_src),1], [ni(i,1,i_src),nj(i,1,i_src),inf] );
                                
                                % Interpolation weights
                                ii = squeeze( ntrp_info{i_src}{1}{1}(1,:,:) ) - i0(i,1,i_src)+1;
                                jj = squeeze( ntrp_info{i_src}{2}{1}(1,:,:) ) - j0(i,1,i_src)+1;
                                ww = squeeze( ntrp_info{i_src}{3}{1}(1,:,:,:) );
                                
                                % New data
                                zeta_out = zeros(size(output_grid.lon_rho,2),size(zeta_in,3));
                                
                                % Interpolate to new grid
                                for j=1:size(zeta_out,1)
                                if(output_grid.mask_rho(1,j)==1)
                                   
                                  % The 4 surrounding data point inidices
                                  iii = ii(j,:);
                                  jjj = jj(j,:);
                                  www = squeeze(ww(j,:,:));
                                  
%                                   % Is this right?
%                                   xxx = input_grid(i_src).lon_rho(iii+i0(i,1,i_src)-1,jjj+j0(i,1,i_src)-1);
%                                   yyy = input_grid(i_src).lat_rho(iii+i0(i,1,i_src)-1,jjj+j0(i,1,i_src)-1);
%                                   plot(xxx,yyy,'.k',output_grid.lon_rho(1,j),output_grid.lat_rho(1,j),'or');
                                  
                                  % The 4 surrounding points
                                  zzz = zeta_in(iii,jjj,:);
                                  
                                  % Do interpolations
                                  zzz = zzz.*repmat(www,[1 1 size(zeta_in,3)]);
                                  zzz = sum(sum(zzz,1),2);
                                  
                                  % Save
                                  zeta_out(j,:) = squeeze(zzz);
                                  
                                  % Clean-up
                                  clear iii jjj www zzz;
                                    
                                end
                                end
                                clear j;
                                
                                % Save to file
                                ncwrite(bry_fileopts.FileName,['zeta_' str_wsen{i}],zeta_out,[1 output_it0]);
                                
                            % . . . . . . . . . . . . . . . . . . . . . . %   
                            case 2  % South
                                
                                zeta_in  = ncread(input_filename, 'zeta', [i0(i,1,i_src),j0(i,1,i_src),1], [ni(i,1,i_src),nj(i,1,i_src),inf] );
                                ii       = squeeze( ntrp_info{i_src}{1}{1}(:,1,:) ) - i0(i,1,i_src)+1;
                                jj       = squeeze( ntrp_info{i_src}{2}{1}(:,1,:) ) - j0(i,1,i_src)+1;
                                ww       = squeeze( ntrp_info{i_src}{3}{1}(:,1,:,:) );
                                zeta_out = zeros(size(output_grid.lon_rho,1),size(zeta_in,3));
                                for j=1:size(zeta_out,1)
                                if(output_grid.mask_rho(j,1)==1)
                                  zeta_out(j,:) = squeeze(sum(sum( zeta_in(ii(j,:),jj(j,:),:) .* repmat(squeeze(ww(j,:,:)),[1 1 size(zeta_in,3)]) ,1),2));
                                end
                                end
                                ncwrite(bry_fileopts.FileName,['zeta_' str_wsen{i}],zeta_out,[1 output_it0]);
                                
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            % East
                            case 3
                                
                                zeta_in  = ncread(input_filename, 'zeta', [i0(i,1,i_src),j0(i,1,i_src),1], [ni(i,1,i_src),nj(i,1,i_src),inf] );
                                ii       = squeeze( ntrp_info{i_src}{1}{1}(end,:,:) ) - i0(i,1,i_src)+1;
                                jj       = squeeze( ntrp_info{i_src}{2}{1}(end,:,:) ) - j0(i,1,i_src)+1;
                                ww       = squeeze( ntrp_info{i_src}{3}{1}(end,:,:,:) );
                                zeta_out = zeros(size(output_grid.lon_rho,2),size(zeta_in,3));
                                for j=1:size(zeta_out,1)
                                if(output_grid.mask_rho(end,j)==1)
                                  zeta_out(j,:) = squeeze(sum(sum( zeta_in(ii(j,:),jj(j,:),:) .* repmat(squeeze(ww(j,:,:)),[1 1 size(zeta_in,3)]) ,1),2));
                                end
                                end
                                ncwrite(bry_fileopts.FileName,['zeta_' str_wsen{i}],zeta_out,[1 output_it0]);
                            
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            % North
                            case 4
                                % NOT YET IMPLEMENTED
                        end
                        
                    end
                    end
                    clear i;
                    
                %% ----------------------------------------------------- %%
                case 'ubar'
                    
                    for i=1:4
                    if(bry_fileopts.WSEN(i)==1)
                        
                        % Second input file
                        input_filenam2 = [src_info{i_src,2} 'vbar/' src_info{i_src,1} '_vbar_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
        
                        % Gather data surrounding boundary
                        switch i
                            
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            case 1  % West
                                
                                % Input data
                                u_in  = ncread(input_filename, 'ubar', [i0(i,2,i_src),j0(i,2,i_src),1], [ni(i,2,i_src),nj(i,2,i_src),inf] );
                                v_in  = ncread(input_filenam2, 'vbar', [i0(i,3,i_src),j0(i,3,i_src),1], [ni(i,3,i_src),nj(i,3,i_src),inf] );
                                
                                % Interp points
                                u_ii = squeeze( ntrp_info{i_src}{1}{2}(1,:,:) ); % - i0(i,2,i_src)+1;
                                u_jj = squeeze( ntrp_info{i_src}{2}{2}(1,:,:) ); % - j0(i,2,i_src)+1;
                                u_ww = squeeze( ntrp_info{i_src}{3}{2}(1,:,:,:) );
                                v_ii = squeeze( ntrp_info{i_src}{1}{3}(1,:,:) ); % - i0(i,3,i_src)+1;
                                v_jj = squeeze( ntrp_info{i_src}{2}{3}(1,:,:) ); %% - j0(i,3,i_src)+1;
                                v_ww = squeeze( ntrp_info{i_src}{3}{3}(1,:,:,:) );
                                
                                % Output variables
                                u_out = zeros(size(output_grid.lon_u,2),size(u_in,3));
                                v_out = zeros(size(output_grid.lon_v,2),size(v_in,3));
                                
                                % Loop through times
                                for j=1:size(u_in,3)
                                    
                                    % For u  
                                    for k=1:size(u_out,1)
                                    
                                        % U
                                        uu_in = u_in(u_ii(k,:)-i0(i,2,i_src)+1,:,j);
                                    
                                        % Grid info
                                        uu_xx = input_grid(i_src).lon_u( u_ii(k,:),u_jj(k,:));
                                        uu_yy = input_grid(i_src).lat_u( u_ii(k,:),u_jj(k,:));
                                        uu_mm = input_grid(i_src).mask_u(u_ii(k,:),u_jj(k,:));
                                        vv_xx = input_grid(i_src).lon_v( v_ii(k,:),v_jj(k,:));
                                        vv_yy = input_grid(i_src).lat_v( v_ii(k,:),v_jj(k,:));
                                        vv_mm = input_grid(i_src).mask_v(v_ii(k,:),v_jj(k,:));
                                    
                                        % Interpolants
                                        uu_ntrplnt = scatteredInterpolant( uu_xx(uu_mm==1), uu_yy(uu_mm==1), uu_in(uu_mm==1), 'linear' );
                                        vv_ntrplnt = scatteredInterpolant( vv_xx(vv_mm==1), vv_yy(uu_mm==1), vv_in(vv_mm==1), 'linear' );
                                    
                                    % Interpolate to each other's grids
                                    uu_onv = uu_ntrplnt(vv_xx,vv_yy);
                                    vv_onu = vv_ntrplnt(uu_xx,uu_yy);
                                    
                                    % Rotate velocities to be eastward, northward
                                    uu_theta = -input_grid(i_src).angle_u(u_ii,u_jj);
                                    vv_theta = -input_grid(i_src).angle_v(v_ii,v_jj);
                                    zz_onu   = (uu_in  + sqrt(-1).*vv_onu) .* exp(sqrt(-1).*uu_theta);
                                    zz_onv   = (uu_onv + sqrt(-1).*vv_in ) .* exp(sqrt(-1).*vv_theta);
                                    
                                    % Interpolate to new grid points
                                    for k=1:size(u_out,1)
                                    if(output_grid.mask_u(1,k)==1)
                                        zzz        = zz_onu( u_ii-i0(i,2,i_src)+1 , u_jj-j0(i,2,i_src)+1 );
                                        www        = squeeze(u_ww(k,:,:));
                                        u_out(k,j) = real( sum(sum( zzz.*www )) * exp(sqrt(-1).*ouput_grid.angle_u(1,k)) );
                                    end
                                    end
                                    for k=1:size(v_out,1)
                                    if(output_grid.mask_v(1,k)==1)
                                        zzz        = zz_onv( v_ii-i0(i,3,i_src)+1 , v_jj-j0(i,3,i_src)+1 );
                                        www        = squeeze(v_ww(k,:,:));
                                        v_out(k,j) = real( sum(sum( zzz.*www )) * exp(sqrt(-1).*ouput_grid.angle_v(1,k)) );
                                    end
                                    end
                                    
                                    % Clean-up
                                    clear k uu_* vv_* zz_* zzz www;
                                    
                                end
                                
                                % Save
                                ncwrite(bry_fileopts.FileName,['u_' str_wsen{i}],u_out,[1 output_it0]);
                                ncwrite(bry_fileopts.FileName,['v_' str_wsen{i}],v_out,[1 output_it0]);
                                
                                % Clean-up
                                clear u_* v_* j
                            
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            case 2  % South
                                
                                u_in  = ncread(input_filename, 'ubar', [i0(i,2,i_src),j0(i,2,i_src),1], [ni(i,2,i_src),nj(i,2,i_src),inf] );
                                v_in  = ncread(input_filename, 'vbar', [i0(i,3,i_src),j0(i,3,i_src),1], [ni(i,3,i_src),nj(i,3,i_src),inf] );
                                u_ii = squeeze( ntrp_info{i_src}{1}{2}(:,1,:) ); % - i0(i,2,i_src)+1;
                                u_jj = squeeze( ntrp_info{i_src}{2}{2}(:,1,:) ); % - j0(i,2,i_src)+1;
                                u_ww = squeeze( ntrp_info{i_src}{3}{2}(:,1,:,:) );
                                v_ii = squeeze( ntrp_info{i_src}{1}{3}(:,1,:) ); % - i0(i,3,i_src)+1;
                                v_jj = squeeze( ntrp_info{i_src}{2}{3}(:,1,:) ); %% - j0(i,3,i_src)+1;
                                v_ww = squeeze( ntrp_info{i_src}{3}{3}(:,1,:,:) );
                                u_out = zeros(size(output_grid.lon_u,1),size(u_in,3));
                                v_out = zeros(size(output_grid.lon_v,1),size(v_in,3));
                                for j=1:size(u_in,3)
                                    uu_in = u_in(:,:,j);
                                    vv_in = v_in(:,:,j);
                                    uu_xx = input_grid(i_src).lon_u( u_ii,u_jj);
                                    uu_yy = input_grid(i_src).lat_u( u_ii,u_jj);
                                    uu_mm = input_grid(i_src).mask_u(u_ii,u_jj);
                                    vv_xx = input_grid(i_src).lon_v( v_ii,v_jj);
                                    vv_yy = input_grid(i_src).lat_v( v_ii,v_jj);
                                    vv_mm = input_grid(i_src).mask_v(v_ii,v_jj);
                                    uu_ntrplnt = scatteredInterpolant( uu_xx(uu_mm==1), uu_yy(uu_mm==1), uu_in(uu_mm==1), 'linear' );
                                    vv_ntrplnt = scatteredInterpolant( vv_xx(vv_mm==1), vv_yy(uu_mm==1), vv_in(vv_mm==1), 'linear' );
                                    uu_onv = uu_ntrplnt(vv_xx,vv_yy);
                                    vv_onu = vv_ntrplnt(uu_xx,uu_yy);
                                    uu_theta = -input_grid(i_src).angle_u(u_ii,u_jj);
                                    vv_theta = -input_grid(i_src).angle_v(v_ii,v_jj);
                                    zz_onu   = (uu_in  + sqrt(-1).*vv_onu) .* exp(sqrt(-1).*uu_theta);
                                    zz_onv   = (uu_onv + sqrt(-1).*vv_in ) .* exp(sqrt(-1).*vv_theta);
                                    for k=1:size(u_out,1)
                                    if(output_grid.mask_u(k,1)==1)
                                        zzz        = zz_onu( u_ii-i0(i,2,i_src)+1 , u_jj-j0(i,2,i_src)+1 );
                                        www        = squeeze(u_ww(k,:,:));
                                        u_out(k,j) = real( sum(sum( zzz.*www )) * exp(sqrt(-1).*ouput_grid.angle_u(k,1)) );
                                    end
                                    end
                                    for k=1:size(v_out,1)
                                    if(output_grid.mask_v(k,1)==1)
                                        zzz        = zz_onv( v_ii-i0(i,3,i_src)+1 , v_jj-j0(i,3,i_src)+1 );
                                        www        = squeeze(v_ww(k,:,:));
                                        v_out(k,j) = real( sum(sum( zzz.*www )) * exp(sqrt(-1).*ouput_grid.angle_v(k,1)) );
                                    end
                                    end
                                end
                                ncwrite(bry_fileopts.FileName,['u_' str_wsen{i}],u_out,[1 output_it0]);
                                ncwrite(bry_fileopts.FileName,['v_' str_wsen{i}],v_out,[1 output_it0]);
                                
                            % . . . . . . . . . . . . . . . . . . . . . . %    
                            case 3  % East
                                
                                % Input data
                                u_in  = ncread(input_filename, 'ubar', [i0(i,2,i_src),j0(i,2,i_src),1], [ni(i,2,i_src),nj(i,2,i_src),inf] );
                                v_in  = ncread(input_filename, 'vbar', [i0(i,3,i_src),j0(i,3,i_src),1], [ni(i,3,i_src),nj(i,3,i_src),inf] );
                                u_ii = squeeze( ntrp_info{i_src}{1}{2}(end,:,:) ); % - i0(i,2,i_src)+1;
                                u_jj = squeeze( ntrp_info{i_src}{2}{2}(end,:,:) ); % - j0(i,2,i_src)+1;
                                u_ww = squeeze( ntrp_info{i_src}{3}{2}(end,:,:,:) );
                                v_ii = squeeze( ntrp_info{i_src}{1}{3}(end,:,:) ); % - i0(i,3,i_src)+1;
                                v_jj = squeeze( ntrp_info{i_src}{2}{3}(end,:,:) ); %% - j0(i,3,i_src)+1;
                                v_ww = squeeze( ntrp_info{i_src}{3}{3}(end,:,:,:) );
                                u_out = zeros(size(output_grid.lon_u,2),size(u_in,3));
                                v_out = zeros(size(output_grid.lon_v,2),size(v_in,3));
                                for j=1:size(u_in,3)
                                    uu_in = u_in(:,:,j);
                                    vv_in = v_in(:,:,j);
                                    uu_xx = input_grid(i_src).lon_u( u_ii,u_jj);
                                    uu_yy = input_grid(i_src).lat_u( u_ii,u_jj);
                                    uu_mm = input_grid(i_src).mask_u(u_ii,u_jj);
                                    vv_xx = input_grid(i_src).lon_v( v_ii,v_jj);
                                    vv_yy = input_grid(i_src).lat_v( v_ii,v_jj);
                                    vv_mm = input_grid(i_src).mask_v(v_ii,v_jj);
                                    uu_ntrplnt = scatteredInterpolant( uu_xx(uu_mm==1), uu_yy(uu_mm==1), uu_in(uu_mm==1), 'linear' );
                                    vv_ntrplnt = scatteredInterpolant( vv_xx(vv_mm==1), vv_yy(uu_mm==1), vv_in(vv_mm==1), 'linear' );
                                    uu_onv = uu_ntrplnt(vv_xx,vv_yy);
                                    vv_onu = vv_ntrplnt(uu_xx,uu_yy);
                                    uu_theta = -input_grid(i_src).angle_u(u_ii,u_jj);
                                    vv_theta = -input_grid(i_src).angle_v(v_ii,v_jj);
                                    zz_onu   = (uu_in  + sqrt(-1).*vv_onu) .* exp(sqrt(-1).*uu_theta);
                                    zz_onv   = (uu_onv + sqrt(-1).*vv_in ) .* exp(sqrt(-1).*vv_theta);
                                    for k=1:size(u_out,1)
                                    if(output_grid.mask_u(end,k)==1)
                                        zzz        = zz_onu( u_ii-i0(i,2,i_src)+1 , u_jj-j0(i,2,i_src)+1 );
                                        www        = squeeze(u_ww(k,:,:));
                                        u_out(k,j) = real( sum(sum( zzz.*www )) * exp(sqrt(-1).*ouput_grid.angle_u(end,k)) );
                                    end
                                    end
                                    for k=1:size(v_out,1)
                                    if(output_grid.mask_v(end,k)==1)
                                        zzz        = zz_onv( v_ii-i0(i,3,i_src)+1 , v_jj-j0(i,3,i_src)+1 );
                                        www        = squeeze(v_ww(k,:,:));
                                        v_out(k,j) = real( sum(sum( zzz.*www )) * exp(sqrt(-1).*ouput_grid.angle_v(end,k)) );
                                    end
                                    end
                                end
                                ncwrite(bry_fileopts.FileName,['u_' str_wsen{i}],u_out,[1 output_it0]);
                                ncwrite(bry_fileopts.FileName,['v_' str_wsen{i}],v_out,[1 output_it0]);
                                
                            % . . . . . . . . . . . . . . . . . . . . . . %    
                            case 4  % North
                        end
                        
                    end
                    end
                    clear i;

                %% ----------------------------------------------------- %%    
                case {'temp','salt'}
                    
                    for i=1:4
                    if(bry_fileopts.WSEN(i)==1)
                        
                        % Gather data surrounding boundary
                        switch i
                            
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            case 1  % West
                                
                                % Input data
                                data_in  = ncread(input_filename, var_to_get{iv,1}, [i0(i,1,i_src),j0(i,1,i_src),1,1], [ni(i,1,i_src),nj(i,1,i_src),inf,inf] );
                                
                                % Interpolation weights
                                ii = squeeze( ntrp_info{i_src}{1}{1}(1,:,:) ) - i0(i,1,i_src)+1;
                                jj = squeeze( ntrp_info{i_src}{2}{1}(1,:,:) ) - j0(i,1,i_src)+1;
                                ww = squeeze( ntrp_info{i_src}{3}{1}(1,:,:,:) );
                                
                                % New data
                                data_out = zeros(size(output_grid.lon_rho,2),numel(output_grid.s_rho),size(data_in,4));
                                
                                % Interpolate to new grid
                                for j=1:size(data_out,1)
                                if(output_grid.mask_rho(1,j)==1)
                                   
                                  % The 4 surrounding data point inidices
                                  iii = ii(j,:);
                                  jjj = jj(j,:);
                                  www = squeeze(ww(j,:,:));
                                  
%                                   % Is this right? (YES!)
%                                   xxx = input_grid(i_src).lon_rho(iii+i0(i,1,i_src)-1,jjj+j0(i,1,i_src)-1);
%                                   yyy = input_grid(i_src).lat_rho(iii+i0(i,1,i_src)-1,jjj+j0(i,1,i_src)-1);
%                                   plot(xxx,yyy,'.k',output_grid.lon_rho(1,j),output_grid.lat_rho(1,j),'or');
                                  
                                  % The 4 surrounding points
                                  zzz = zeta_in(iii,jjj,:,:);
                                  
                                  % Interpolate data in depth
                                  znew = zeros(size(zzz,1),size(zzz,2),numel(output_grid.s_rho),size(zzz,4));
                                  for a=1:2
                                  for b=1:2
                                      znew(a,b,:,:) = interp1( input_grid(i_src).s_rho, squeeze(zzz(a,b,:,:)), output_grid.s_rho, 'linear' );
                                  end
                                  end
                                  
                                  % Do interpolations
                                  zzz = zzz.*repmat(www,[1 1 size(znew,3) size(znew,4)]);
                                  zzz = sum(sum(zzz,1),2);
                                  
                                  % Save
                                  data_out(j,:,:) = squeeze(zzz);
                                  
                                  % Clean-up
                                  clear iii jjj www zzz;
                                    
                                end
                                end
                                clear j;
                                
                                % Save to file
                                ncwrite(bry_fileopts.FileName,[ var_to_get{iv,1} '_' str_wsen{i}],data_out,[1 1 output_it0]);
                                
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            case 2  % South
                               
                                % Input data
                                data_in  = ncread(input_filename, var_to_get{iv,1}, [i0(i,1,i_src),j0(i,1,i_src),1,1], [ni(i,1,i_src),nj(i,1,i_src),inf,inf] );
                                ii       = squeeze( ntrp_info{i_src}{1}{1}(:,1,:) ) - i0(i,1,i_src)+1;
                                jj       = squeeze( ntrp_info{i_src}{2}{1}(:,1,:) ) - j0(i,1,i_src)+1;
                                ww       = squeeze( ntrp_info{i_src}{3}{1}(:,1,:,:) );
                                data_out = zeros(size(output_grid.lon_rho,1),numel(output_grid.s_rho),size(data_in,4));
                                for j=1:size(data_out,1);   if(output_grid.mask_rho(1,j)==1)
                                    zzz = zeta_in(ii(j,:),jj(j,:),:,:);
                                    znew = zeros(size(zzz,1),size(zzz,2),numel(output_grid.s_rho),size(zzz,4));
                                    for a=1:2;  for b=1:2
                                        znew(a,b,:,:) = interp1( input_grid(i_src).s_rho, squeeze(zzz(a,b,:,:)), output_grid.s_rho, 'linear' );
                                    end;        end
                                    zzz = zzz.*repmat(squeeze(ww(j,:,:)),[1 1 size(znew,3) size(znew,4)]);
                                    data_out(j,:,:) = squeeze(sum(sum(zzz,1),2));
                                end;                        end
                                ncwrite(bry_fileopts.FileName,[ var_to_get{iv,1} '_' str_wsen{i}],data_out,[1 1 output_it0]);
                                
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            case 3  % East
                                
                                % Input data
                                data_in  = ncread(input_filename, var_to_get{iv,1}, [i0(i,1,i_src),j0(i,1,i_src),1,1], [ni(i,1,i_src),nj(i,1,i_src),inf,inf] );
                                ii       = squeeze( ntrp_info{i_src}{1}{1}(end,:,:) ) - i0(i,1,i_src)+1;
                                jj       = squeeze( ntrp_info{i_src}{2}{1}(end,:,:) ) - j0(i,1,i_src)+1;
                                ww       = squeeze( ntrp_info{i_src}{3}{1}(end,:,:,:) );
                                data_out = zeros(size(output_grid.lon_rho,2),numel(output_grid.s_rho),size(data_in,4));
                                for j=1:size(data_out,1);   if(output_grid.mask_rho(1,j)==1)
                                    zzz = zeta_in(ii(j,:),jj(j,:),:,:);
                                    znew = zeros(size(zzz,1),size(zzz,2),numel(output_grid.s_rho),size(zzz,4));
                                    for a=1:2;  for b=1:2
                                        znew(a,b,:,:) = interp1( input_grid(i_src).s_rho, squeeze(zzz(a,b,:,:)), output_grid.s_rho, 'linear' );
                                    end;        end
                                    zzz = zzz.*repmat(squeeze(ww(j,:,:)),[1 1 size(znew,3) size(znew,4)]);
                                    data_out(j,:,:) = squeeze(sum(sum(zzz,1),2));
                                end;                        end
                                ncwrite(bry_fileopts.FileName,[ var_to_get{iv,1} '_' str_wsen{i}],data_out,[1 1 output_it0]);
                                
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            case 4  % North
                        end
                    end
                    end
                    clear i;
                        
                %% ----------------------------------------------------- %%   
                case 'u'
           
                   % Second input file
                   input_filenam2 = [src_info{i_src,2} 'v/' src_info{i_src,1} '_v_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
        
                    
                   for i=1:4
                    if(bry_fileopts.WSEN(i)==1)
                        
                        % Gather data surrounding boundary
                        switch i
                            
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            case 1  % West
                                
                                % Input data
                                u_in  = ncread(input_filename, 'u', [i0(i,2,i_src),j0(i,2,i_src),1,1], [ni(i,2,i_src),nj(i,2,i_src),inf,inf] );
                                v_in  = ncread(input_filenam2, 'v', [i0(i,3,i_src),j0(i,3,i_src),1,1], [ni(i,3,i_src),nj(i,3,i_src),inf,inf] );
                                
                                % Interp points
                                u_ii = squeeze( ntrp_info{i_src}{1}{2}(1,:,:) ); % - i0(i,2,i_src)+1;
                                u_jj = squeeze( ntrp_info{i_src}{2}{2}(1,:,:) ); % - j0(i,2,i_src)+1;
                                u_ww = squeeze( ntrp_info{i_src}{3}{2}(1,:,:,:) );
                                v_ii = squeeze( ntrp_info{i_src}{1}{3}(1,:,:) ); % - i0(i,3,i_src)+1;
                                v_jj = squeeze( ntrp_info{i_src}{2}{3}(1,:,:) ); %% - j0(i,3,i_src)+1;
                                v_ww = squeeze( ntrp_info{i_src}{3}{3}(1,:,:,:) );
                                
                                % Output variables
                                u_out = zeros(size(output_grid.lon_u,2),numel(output_grid.s_rho),size(u_in,4));
                                v_out = zeros(size(output_grid.lon_v,2),numel(output_grid.s_rho),size(v_in,4));
                                
                                % Loop through times
                                for j=1:size(u_in,4)
                                    
                                    % U's and V's
                                    uu_in = u_in(:,:,:,j);
                                    vv_in = v_in(:,:,:,j);
                                    
                                    % Grid info
                                    uu_xx = input_grid(i_src).lon_u( u_ii,u_jj);
                                    uu_yy = input_grid(i_src).lat_u( u_ii,u_jj);
                                    uu_mm = input_grid(i_src).mask_u(u_ii,u_jj);
                                    vv_xx = input_grid(i_src).lon_v( v_ii,v_jj);
                                    vv_yy = input_grid(i_src).lat_v( v_ii,v_jj);
                                    vv_mm = input_grid(i_src).mask_v(v_ii,v_jj);
                                    
                                    % Interpolate to each other's grid at each depth
                                    uu_onv = zeros(size(vv_in));
                                    vv_onu = zeros(size(uu,in));
                                    for k=1:size(uu_in,3)
                                        
                                        % Datums
                                        uuu = uu_in(:,:,k);
                                        vvv = vv_in(:,:,k);
                                        
                                        % Interpolants
                                        uu_ntrplnt = scatteredInterpolant( uu_xx(uu_mm==1), uu_yy(uu_mm==1), uuu(uu_mm==1), 'linear' );
                                        vv_ntrplnt = scatteredInterpolant( vv_xx(vv_mm==1), vv_yy(uu_mm==1), vvv(vv_mm==1), 'linear' );
                                    
                                        % Interpolate to each other's grids
                                        uu_onv(:,:,k) = uu_ntrplnt(vv_xx,vv_yy);
                                        vv_onu(:,:,k) = vv_ntrplnt(uu_xx,uu_yy);
                                        
                                    end
                                    
                                    % Rotate velocities to be eastward, northward
                                    uu_theta = -repmat( input_grid(i_src).angle_u(u_ii,u_jj) , [1 1 size(uu_in,3)]);
                                    vv_theta = -repmat( input_grid(i_src).angle_v(v_ii,v_jj) , [1 1 size(vv_in,3)]);
                                    zz_onu   = (uu_in  + sqrt(-1).*vv_onu) .* exp(sqrt(-1).*uu_theta);
                                    zz_onv   = (uu_onv + sqrt(-1).*vv_in ) .* exp(sqrt(-1).*vv_theta);
                                    
                                    % Interpolate in depth
                                    znew_onu = zeros(size(zz_onu,1),size(zz_onu,2),numel(output_grid.s_rho));
                                    znew_onv = zeros(size(zz_onv,1),size(zz_onv,2),numel(output_grid.s_rho));
                                    for a=1:size(zz_onu,1)
                                    for b=1:size(zz_onu,2)
                                        znew_onu(a,b,:) = interp1( input_grid(i_src).s_rho, squeeze(zz_onu(a,b,:)), output_grid.s_rho, 'linear' );
                                    end
                                    end
                                    for a=1:size(zz_onv,1)
                                    for b=1:size(zz_onv,2)
                                        znew_onv(a,b,:) = interp1( input_grid(i_src).s_rho, squeeze(zz_onv(a,b,:)), output_grid.s_rho, 'linear' );
                                    end
                                    end
                                    
                                    % Interpolate to new grid points
                                    for k=1:size(u_out,1)
                                    if(output_grid.mask_u(1,k)==1)
                                        for m=1:size(u_out,3)
                                            zzz          = znew_onu( u_ii-i0(i,2,i_src)+1 , u_jj-j0(i,2,i_src)+1 , m );
                                            www          = squeeze(u_ww(k,:,:));
                                            u_out(k,m,j) = real( sum(sum( zzz.*www )) * exp(sqrt(-1).*ouput_grid.angle_u(1,k)) );
                                        end
                                    end
                                    end
                                    for k=1:size(v_out,1)
                                    if(output_grid.mask_v(1,k)==1)
                                        for m=1:size(v_out,3)
                                            zzz          = znew_onv( v_ii-i0(i,3,i_src)+1 , v_jj-j0(i,3,i_src)+1 , m );
                                            www          = squeeze(v_ww(k,:,:));
                                            v_out(k,m,j) = real( sum(sum( zzz.*www )) * exp(sqrt(-1).*ouput_grid.angle_v(1,k)) );
                                        end
                                    end
                                    end
                                    
                                    % Clean-up
                                    clear k uu_* vv_* zz_* zzz www;
                                    
                                end
                                
                                % Save
                                ncwrite(bry_fileopts.FileName,['u_' str_wsen{i}],u_out,[1 output_it0]);
                                ncwrite(bry_fileopts.FileName,['v_' str_wsen{i}],v_out,[1 output_it0]);
                                
                                % Clean-up
                                clear u_* v_* j
                            
                            % . . . . . . . . . . . . . . . . . . . . . . %
                            case 2  % South
                                pause(1);
                                
                            case 3
                                
                                
                            case 4
                                
                        end
                        
                    end
                   end
                    
                otherwise
                    %NOTHING
            end

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