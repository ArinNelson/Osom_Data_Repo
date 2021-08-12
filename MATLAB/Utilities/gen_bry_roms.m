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
Switch(1) = 1;      % Compute interpolation info things
Switch(2) = 1;      % Init the boundary file
Switch(3) = 1;      % Gather & interpolate data to boundaries

% Time Options
date_start = [2018,01,01];	% Start year, month, day to gather data for
date_end   = [2018,12,31];	% End   year, month, day to gather data for
date_plus  = false;         % Have last entry be the first timestep of the day after date_end

% Grid options
output_gridfile = 'C:/Library/ROMS_Stuff/Resources/ngbay_grd.nc';                                   % ROMS grid file
output_infofile = 'C:/Library/ROMS_Stuff/Resources/forcefiles_riroms/riroms_bry_2006_dt06hrs.nc';   % Contains vertical grid info

% Other options
smooth_factor = 0;      % Low-pass filter with this window (in days).   e.g., to de-tide NECOFS, Dave Ullman uses 36hrs
bry_fileopts = struct;	% Options for save file
bry_fileopts.FileName = 'OSOM_bry_2018_DOPPIO_best.nc';
bry_fileopts.Title    = 'BOUNDARY forcing file: DOPPIO (best) source';
bry_fileopts.WSEN     = [1 1 1 0];      % Include west, south, east, and/or north boundaries

% Variable info (Var Label , Source , dimensions (2 or 3), grid type 
var_to_get = {'zeta',  'DOPPIO_best',   2,  1;  ...
              'ubar',  'DOPPIO_best',   2,  2;  ...
              'vbar',  'DOPPIO_best',   2,  3;  ...
              'u',     'DOPPIO_best',   3,  2;  ...
              'v',     'DOPPIO_best',   3,  3;  ...
              'temp',  'DOPPIO_best',   3,  1;  ...
              'salt',  'DOPPIO_best',   3,  1;  ...
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
        ntrp_i = cell(3,1);
        ntrp_j = cell(3,1);
        ntrp_k = cell(3,1);
        
        % For rho
        [xx,yy] = grn2eqa(input_grid(is).lat_rho,input_grid(is).lon_rho);
        mm      = input_grid(is).mask_rho;
        [xq,yq] = grn2eqa(output_grid.lat_rho,output_grid.lon_rho);
        mq      = output_grid.mask_rho;
        [ntrp_i{1}, ntrp_j{1}, ntrp_w{1}] = grid_interpolant(xx,yy,xq,yq,mm,mq);
        
        % For u
        [xx,yy] = grn2eqa(input_grid(is).lat_u,  input_grid(is).lon_u  );
        mm      = input_grid(is).mask_u;
        [xq,yq] = grn2eqa(output_grid.lat_u  ,output_grid.lon_u  );
        mq      = output_grid.mask_u;
        [ntrp_i{2}, ntrp_j{2}, ntrp_w{2}] = grid_interpolant(xx,yy,xq,yq,mm,mq);

        % For v
        [xx,yy] = grn2eqa(input_grid(is).lat_v,  input_grid(is).lon_v  );
        mm      = input_grid(is).mask_v;
        [xq,yq] = grn2eqa(output_grid.lat_v  ,output_grid.lon_v  );
        mq      = output_grid.mask_v;
        [ntrp_i{3}, ntrp_j{3}, ntrp_w{3}] = grid_interpolant(xx,yy,xq,yq,mm,mq);
        
        % Save to variable
        interp_info = {ntrp_i, ntrp_j, ntrp_w};
        ntrp_info{is} = interp_info;
         
        % Also backup to file
        save(['bry_ntrp_info_' input_src{is} '.mat'],'interp_info');
        clear interp_info;
         
    else
        
        % Load previously-defined interpolation info
        tmp = load(['bry_ntrp_info_' input_src{is} '.mat']);  
        ntrp_info{is} = tmp.interp_info;
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
        test     = ceil( test(end) + datenum(date_start(1),date_start(2),date_start(3)) );
        year_on  = year(test);
        month_on = month(test);
    end
    clear test;

    % WSEN strings
    str_wsen = {'west','south','east','north'};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Loop through years & months
    while( datenum(year_on,month_on,1) <= datenum(date_end(1),date_end(2),1) )
    for i=1:4   % WSEN
        
    %=====================================================================%    
    if(bry_fileopts.WSEN(i)==1)
    disp(['Gathering data for ' num2str(year_on) '/' sprintf('%0.2d',month_on) ', boundary ' str_wsen{i} '...']);
    
        % Output grid
        outgrid	= grid_border(output_grid,i);
        
        % Loop through variables
        for iv=1:n_var
            
            %.............................................................%
            
            % Input source
            i_src = find( strcmp( src_info(:,1), var_to_get{iv,2} )==1 );
            
            % Input grid
            ingrid = input_grid(i_src);
            
            % Grid variables & interpolation info
            switch var_to_get{iv,4}
                case 1
                    out_mask  = outgrid.mask_rho;   in_mask  = ingrid.mask_rho;
                    ntrpinfo = ntrp_info{i_src};
                    for ii=1:3
                        ntrpinfo{ii} = ntrpinfo{ii}{1};
                        switch i
                            case 1;     ntrpinfo{ii} = squeeze( ntrpinfo{ii}(1,:,:,:)   );
                            case 2;     ntrpinfo{ii} = squeeze( ntrpinfo{ii}(:,1,:,:)   );
                            case 3;     ntrpinfo{ii} = squeeze( ntrpinfo{ii}(end,:,:,:) );
                            case 4;     ntrpinfo{ii} = squeeze( ntrpinfo{ii}(:,end,:,:) );
                        end
                    end
                otherwise
                    out_lon_u   = outgrid.lon_u;    in_lon_u   = ingrid.lon_u;
                    out_lon_v   = outgrid.lon_v;    in_lon_v   = ingrid.lon_v;
                    out_lat_u   = outgrid.lat_u;    in_lat_u   = ingrid.lat_u;
                    out_lat_v   = outgrid.lat_v;    in_lat_v   = ingrid.lat_v;
                    out_mask_u  = outgrid.mask_u;   in_mask_u  = ingrid.mask_u;
                    out_mask_v  = outgrid.mask_v;   in_mask_v  = ingrid.mask_v;
                    out_angle_u = outgrid.angle_u;  in_angle_u = ingrid.angle_u;
                    out_angle_v = outgrid.angle_v;  in_angle_v = ingrid.angle_v;
                    ntrpinfo_u = ntrp_info{i_src};  for j=1:3;  ntrpinfo_u{j} = ntrpinfo_u{j}{2};   end
                    ntrpinfo_v = ntrp_info{i_src};  for j=1:3;  ntrpinfo_v{j} = ntrpinfo_v{j}{3};   end
                    for ii=1:3
                        switch i
                            case 1;     ntrpinfo_u{ii} = squeeze( ntrpinfo_u{ii}(1,:,:,:)   );  ntrpinfo_v{ii} = squeeze( ntrpinfo_v{ii}(1,:,:,:)   );
                            case 2;     ntrpinfo_u{ii} = squeeze( ntrpinfo_u{ii}(:,1,:,:)   );  ntrpinfo_v{ii} = squeeze( ntrpinfo_v{ii}(:,1,:,:)   );
                            case 3;     ntrpinfo_u{ii} = squeeze( ntrpinfo_u{ii}(end,:,:,:) );  ntrpinfo_v{ii} = squeeze( ntrpinfo_v{ii}(end,:,:,:) );
                            case 4;     ntrpinfo_u{ii} = squeeze( ntrpinfo_u{ii}(:,end,:,:) );  ntrpinfo_v{ii} = squeeze( ntrpinfo_v{ii}(:,end,:,:) );
                        end
                    end
            end
            
            % Input file name
            input_file = [src_info{i_src,2} var_to_get{iv,1} '/' src_info{i_src,1} '_' var_to_get{iv,1} '_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
        
            % If first variable and border, get time info
            if(exist('output_it0','var')~=1)
                
                % Input time info
                input_time = ncread(input_file,'time');
                input_it0  = find(input_time >= datenum(date_start(1),date_start(2),date_start(3)),1,'first');
                input_itf  = find(input_time <= datenum(date_end(1),  date_end(2),  date_end(3))+1,1,'last');
                input_nt   = input_itf-input_it0+1;
                
                % Where to start saving to output file
                output_time = ncread(bry_fileopts.FileName,'bry_time');
                output_time(output_time>1e36) = [];
                if(isempty(output_time)==1)
                    output_it0 = 1;
                else
                    output_it0 = numel(output_time)+1;
                end
                clear output_time;
                
            end
            
            % Also get variable sampling info
            switch var_to_get{iv,4}
                case 1
                    ii0  = min(ntrpinfo{1}(ntrpinfo{1}>0));     nii = max(ntrpinfo{1}(:))-ii0+1;        iii = ii0:(ii0+nii-1);
                    jj0  = min(ntrpinfo{2}(ntrpinfo{2}>0));     njj = max(ntrpinfo{2}(:))-jj0+1;        jjj = jj0:(jj0+njj-1);
                    in_mask  = in_mask(iii,jjj);
                otherwise
                    ii0_u = min(ntrpinfo_u{1}(ntrpinfo_u{1}>0));  nii_u = max(ntrpinfo_u{1}(:))-ii0_u+1;	iii_u = ii0_u : (ii0_u+nii_u-1);
                    jj0_u = min(ntrpinfo_u{2}(ntrpinfo_u{2}>0));  njj_u = max(ntrpinfo_u{2}(:))-jj0_u+1;	jjj_u = jj0_u : (jj0_u+njj_u-1);
                    in_lon_u    = in_lon_u(  iii_u,jjj_u);
                    in_lat_u    = in_lat_u(  iii_u,jjj_u);
                    in_mask_u   = in_mask_u( iii_u,jjj_u);
                    in_angle_u  = in_angle_u(iii_u,jjj_u);
                    ii0_v = min(ntrpinfo_v{1}(ntrpinfo_v{1}>0));  nii_v = max(ntrpinfo_v{1}(:))-ii0_v+1;	iii_v = ii0_v : (ii0_v+nii_v-1);
                    jj0_v = min(ntrpinfo_v{2}(ntrpinfo_v{2}>0));  njj_v = max(ntrpinfo_v{2}(:))-jj0_v+1;	jjj_v = jj0_v : (jj0_v+njj_v-1);
                    in_lon_v    = in_lon_v(  iii_v,jjj_v);
                    in_lat_v    = in_lat_v(  iii_v,jjj_v);
                    in_mask_v   = in_mask_v( iii_v,jjj_v);
                    in_angle_v  = in_angle_v(iii_v,jjj_v);
            end
            
            %.............................................................%
            
            switch var_to_get{iv,1}
            
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . %
                case {'zeta','temp','salt'}
                
                    % Load data
                    switch var_to_get{iv,3}
                        case 2
                            data_in     = ncread(input_file,var_to_get{iv,1},[ii0 jj0 1],[nii njj inf]);
                            tmp         = zeros(size(data_in,1),size(data_in,2),1,size(data_in,3));
                            for j=1:size(data_in,3);    tmp(:,:,1,j) = data_in(:,:,j);  end
                            data_in     = tmp;  clear tmp;
                            data_out    = zeros(numel(out_mask),1,input_nt);
                        case 3
                            data_in     = ncread(input_file,var_to_get{iv,1},[ii0 jj0 1 1],[nii njj inf inf]);
                            data_out    = zeros(numel(out_mask),numel(output_grid.s_rho),input_nt);
                    end
                    
                    % Interpolate input data to output grid border
                    for it=1:input_nt
                        
                        % Data at this time step
                        tmp_in = squeeze( data_in(:,:,:,it) );
                        
                        % Interpolate data in depth if needed
                        if(size(tmp_in,3)>1)
                            new_in = zeros(nii,njj,numel(output_grid.s_rho));
                            for a=1:nii
                            for b=1:njj
                            if(in_mask(a,b)==1)
                               new_in(a,b,:) = interp1( input_grid(i_src).s_rho, squeeze(tmp_in(a,b,:)), output_grid.s_rho, 'linear' ); 
                            end
                            end
                            end
                            tmp_in = new_in;    clear new_in;
                        end
                        
                        % Loop through output grid points
                        for j=1:numel(out_mask)
                        if(out_mask(j)==1)
                            
                            % Interpolation info at this point
                            i_ntrp = ntrpinfo{1}(j,:);
                            j_ntrp = ntrpinfo{2}(j,:);
                            w_ntrp = ntrpinfo{3}(j,:,:);
                            
                            % Loop through depths
                            for k=1:size(data_out,2)
                            
                                % Data at this point
                                zz_in = tmp_in(i_ntrp-ii0+1,j_ntrp-jj0+1,k);
                                
                                % Save as weighted sum of the 4 points
                                data_out(j,k,it) = sum( zz_in(:).*w_ntrp(:) ) ./ sum( w_ntrp(:) );
                                
                            end
                            
                            % Clean-up
                            clear i_ntrp j_ntrp w_ntrp k;
                            
                        end
                        end
                        
                        % Clean-up
                        clear tmp_in j;
                        
                    end
                    clear it;
                    
                    % Save to output file
                    switch var_to_get{iv,3}
                        case 2;	ncwrite(bry_fileopts.FileName, [var_to_get{iv,1} '_' str_wsen{i}], squeeze(data_out), [1 output_it0  ]);
                        case 3; ncwrite(bry_fileopts.FileName, [var_to_get{iv,1} '_' str_wsen{i}], data_out,          [1 1 output_it0]);
                    end
                    
                    % Clean-up
                    clear data_in data_out;
                    
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . %
                case {'ubar','u'}
                
                    % Load data
                    input_fil2 = [src_info{i_src,2} var_to_get{iv+1,1} '/' src_info{i_src,1} '_' var_to_get{iv+1,1} '_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
                    switch var_to_get{iv,3}
                        case 2
                            u_onu = ncread(input_file,'ubar',[ii0_u jj0_u 1],[nii_u njj_u inf]);
                            v_onv = ncread(input_fil2,'vbar',[ii0_v jj0_v 1],[nii_v njj_v inf]);
                            tmpu  = zeros(size(u_onu,1),size(u_onu,2),1,size(u_onu,3));
                            tmpv  = zeros(size(v_onv,1),size(v_onv,2),1,size(v_onv,3));
                            for j=1:size(u_onu,3) 
                                tmpu(:,:,1,j) = u_onu(:,:,j);  
                                tmpv(:,:,1,j) = v_onv(:,:,j);
                            end
                            u_onu = tmpu;   v_onv = tmpv;   clear j tmpu tmpv;
                            out_uonu = zeros(numel(out_mask_u),1,input_nt);
                            out_vonv = zeros(numel(out_mask_v),1,input_nt);
                        case 3
                            u_onu = ncread(input_file,'u',[ii0_u jj0_u 1 1],[nii_u njj_u inf inf]);
                            v_onv = ncread(input_fil2,'v',[ii0_v jj0_v 1 1],[nii_v njj_v inf inf]);
                            out_uonu = zeros(numel(out_mask_u),numel(output_grid.s_rho),input_nt);
                            out_vonv = zeros(numel(out_mask_v),numel(output_grid.s_rho),input_nt);
                    end
                    
                    % Interpolate input data to output grid border
                    for it=1:input_nt
                        
                        % Data at this time step
                        uu_onu = squeeze( u_onu(:,:,:,it) );
                        vv_onv = squeeze( v_onv(:,:,:,it) );
                        
                        % Interpolate data in depth if needed
                        if(size(uu_onu,3)>1)
                            
                            % For u
                            new_uu = zeros(size(u_onu,1),size(u_onu,2),numel(output_grid.s_rho));
                            for a=1:size(u_onu,1)
                            for b=1:size(u_onu,2)
                            if(in_mask_u(a,b)==1)
                               new_uu(a,b,:) = interp1( input_grid(i_src).s_rho, squeeze(uu_onu(a,b,:)), output_grid.s_rho, 'linear' ); 
                            end
                            end
                            end
                            uu_onu = new_uu;    clear new_uu;
                            
                            % For v
                            new_vv = zeros(size(v_onv,1),size(v_onv,2),numel(output_grid.s_rho));
                            for a=1:size(v_onv,1)
                            for b=1:size(v_onv,2)
                            if(in_mask_v(a,b)==1)
                               new_vv(a,b,:) = interp1( input_grid(i_src).s_rho, squeeze(vv_onv(a,b,:)), output_grid.s_rho, 'linear' ); 
                            end
                            end
                            end
                            vv_onv = new_vv;    clear new_vv;
                            
                        end
                        
                        % Interpolate to each other's grid at each depth
                        uu_onv = zeros(size(vv_onv));
                        vv_onu = zeros(size(uu_onu));
                        for iz=1:numel(uu_onu,3)
                        
                            % Data
                            tmpu = uu_onu(:,:,iz);
                            tmpv = vv_onv(:,:,iz);
                            
                            % Interpolants
                            uu_ntrplnt = scatteredInterpolant( in_lon_u(in_mask_u==1), in_lat_u(in_mask_u==1), tmpu(in_mask_u==1), 'linear' );
                            vv_ntrplnt = scatteredInterpolant( in_lon_v(in_mask_v==1), in_lat_v(in_mask_v==1), tmpv(in_mask_v==1), 'linear' );
                                    
                            % Interpolate
                            uu_onv(:,:,iz) = uu_ntrplnt(in_lon_v,in_lat_v);
                            vv_onu(:,:,iz) = vv_ntrplnt(in_lon_u,in_lat_u);
                            
                        end
                        clear iz tmpu tmpv;
                        
                        % Rotate velocities to be eastward & northward
                        zz_onu = (uu_onu + sqrt(-1).*vv_onu) .* exp( sqrt(-1) .* -repmat(in_angle_u,[1 1 size(uu_onu,3)]) );
                        zz_onv = (uu_onv + sqrt(-1).*vv_onv) .* exp( sqrt(-1) .* -repmat(in_angle_v,[1 1 size(vv_onv,3)]) );
                        
                        % Interpolate to output grid boundary
                        new_uonu = zeros( size(out_uonu,1), size(out_uonu,2) );
                        new_vonv = zeros( size(out_vonv,1), size(out_vonv,2) );
                        for iz=1:size(out_uonu,2)
                            
                            % For 
                            for iu=1:size(out_uonu,1)   
                            if(ntrpinfo_u{1}(iu)>0 && ntrpinfo_u{2}(iu)>0)   
                                zzz             = zz_onu( ntrpinfo_u{1}(iu)-ii0_u+1, ntrpinfo_u{2}(iu)-jj0_u+1, iz );
                                www             = squeeze(ntrpinfo_u{3}(iu,:,:));
                                new_uonu(iu,iz) = real( sum(sum( zzz.*www )) * exp(sqrt(-1).* out_angle_u(iu)) );
                            end
                            end
                            
                            % For v
                            for iu=1:size(out_vonv,1)
                            if(ntrpinfo_v{1}(iu)>0 && ntrpinfo_v{2}(iu)>0)      
                                zzz             = zz_onv( ntrpinfo_v{1}(iu)-ii0_v+1, ntrpinfo_v{2}(iu)-jj0_v+1, iz );
                                www             = squeeze(ntrpinfo_v{3}(iu,:,:));
                                new_vonv(iu,iz) = imag( sum(sum( zzz.*www )) * exp(sqrt(-1).* out_angle_v(iu)) );
                            end
                            end
                            
                        end
                        clear iz iu zzz www;
                        
                        % Save for this timestep
                        out_uonu(:,:,it) = new_uonu;
                        out_vonv(:,:,it) = new_vonv;
                        
                        % Clean-up
                        clear new_* zz_* uu_* vv_*;
                        
                    end
                    clear it;
                
                    % Save
                    switch var_to_get{iv,3}
                        case 2
                            ncwrite(bry_fileopts.FileName, [var_to_get{iv,  1} '_' str_wsen{i}], squeeze(out_uonu), [1 output_it0  ]);
                            ncwrite(bry_fileopts.FileName, [var_to_get{iv+1,1} '_' str_wsen{i}], squeeze(out_vonv), [1 output_it0  ]);
                        case 3
                            ncwrite(bry_fileopts.FileName, [var_to_get{iv,  1} '_' str_wsen{i}], out_uonu,          [1 1 output_it0  ]);
                            ncwrite(bry_fileopts.FileName, [var_to_get{iv+1,1} '_' str_wsen{i}], out_vonv,          [1 1 output_it0  ]);
                    end
                    
                    % Clean-up
                    clear u_onu v_onv out_uonu out_vonv
                    
            end
                
        end
        clear iv;
        
    %=====================================================================%
    % After all varaibles, save time
    end
	if(i==4)
        output_time = input_time - input_time(1) + datenum(year_on,month_on,1) - datenum(date_start(1),date_start(2),date_start(3));
        ncwrite(bry_fileopts.FileName,'bry_time',output_time,output_it0);
        clear output_time output_it0;
	end

    end
    
    % Onto next year & month
	month_on = month_on + 1;
	if(month_on == 13)
        year_on  = year_on + 1;
        month_on = 1;
	end
    
    end
    clear year_on month_on i;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end