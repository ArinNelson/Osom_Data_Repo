function bry_ntrp_fcn(data_dir, var_name, date_on, ntrp_info, WSEN, save_file, z_old, z_new);           ...

    %---------------------------------------------------------------------%
    % Load data
    
    % If ubar or u, get data file for both u- and v- components
    switch var_name
        case {'ubar','u'}
            
            % Variable names
            u_name    = var_name;
            v_name    = u_name;      
            v_name(1) = 'v';
            
            % Variable files
            u_dir = [data_dir '/' u_name '/*_' num2str(date_on(1)) '_' sprintf('%0.2d',date_on(2)) '.nc'];
            v_dir = [data_dir '/' v_name '/*_' num2str(date_on(1)) '_' sprintf('%0.2d',date_on(2)) '.nc'];
            
            % Look for files
            u_file = ls(u_dir); u_file = [data_dir '/' u_name '/' u_file];
            v_file = ls(v_dir); v_file = [data_dir '/' v_name '/' v_file];
            
            % Read in data
            u_raw = ncread(u_file,u_name);
            v_raw = ncread(v_file,v_name);
            
        case {'zeta','temp','salt'}

            % File to look for
            x_dir = [data_dir '/' var_name '/*_' num2str(date_on(1)) '_' sprintf('%0.2d',date_on(2)) '.nc'];
            
            % File
            x_file = ls(x_dir); z_file = [data_dir '/' var_name '/' x_file];
            
            % Load data
            x_raw = ncread(x_file,var_name);
            
        otherwise
            error(['Unknown or unsupported variable: ' var_name]);
            
    end
    
    %---------------------------------------------------------------------%
    % If 2D data, set so depth dimension is 1
    
    % If velocity variables, need to change 2 variables.  if zeta, just 1
    switch var_name
        
        case 'ubar'
            u_raw_3d = zeros(size(u_raw,1),size(u_raw,2),1,size(u_raw,3));
            v_raw_3d = zeros(size(v_raw,1),size(v_raw,2),1,size(v_raw,3));
            for i=1:size(u_raw,3)
                u_raw_3d(:,:,1,i) = u_raw(:,:,i);
                v_raw_3d(:,:,1,i) = v_raw(:,:,i);
            end
            u_raw = u_raw_3d;   clear u_raw_3d;
            v_raw = v_raw_3d;   clear v_raw_3d;
            
        case 'zeta'
            x_raw_3d = zeros(size(x_raw,1),size(x_raw,2),1,size(x_raw,3));
            for i=1:size(x_raw,3)
                x_raw_3d(:,:,1,i) = x_raw(:,:,i);
            end
            x_raw = x_raw_3d;   clear x_raw_3d;
            
        otherwise
            % Data is already 3D, so do nothing
    end

    %---------------------------------------------------------------------%
    % Interpolate data in depth
    
    % If it's velocity, need to do for 2 variables.
    switch var_name
        
        case {'ubar','u'}
        if(size(u_raw,3)>1)
            
           % New variable
           u_new = zeros(size(u_raw,1),size(u_raw,2),numel(z_new),size(u_raw,4));
           v_new = zeros(size(v_raw,1),size(v_raw,2),numel(z_new),size(v_raw,4));
           for i=1:size(v_raw,1)
           for j=1:size(u_raw,2)
               
             % For u 
             if(i<=size(u_raw,1))
             if(grid(2).mask(i,j)==1)
             for k=1:size(u_raw,4)
                u_new(i,j,:,k) = interp1( z_old, squeeze(u_raw(i,j,:,k)), z_new, 'linear' );
             end
             end
             end

             % for v
             if(j<=size(v_raw,2))
             if(grid(3).mask(i,j)==1)
             for k=1:size(v_raw,4)
                v_new(i,j,:,k) = interp1( z_old, squeeze(v_raw(i,j,:,k)), z_new, 'linear' );
             end
             end
             end
               
           end
           end
           u_raw = u_new;   clear u_new;
           v_raw = v_new;   clear v_new;
            
        end
            
        case {'zeta','temp','salt'}
        if(size(z_raw,3)>1)
            
            % New variable
            x_new = zeros(size(x_raw,1),size(x_raw,2),numel(z_new),size(x_raw,4));
            for i=1:size(x_raw,1)
            for j=1:size(x_raw,2)
            if(grid(1).mask(i,j)==1)
            for k=1:size(x_raw,3)
                x_new(i,j,:,k) = interp1( z_old, squeeze(x_old(i,j,:,k)), z_new, 'linear' );
            end
            end
            end
            end
            x_raw = x_new;  clear x_new;
            
        end
            
        otherwise
            % ?
    end
    
    %---------------------------------------------------------------------%
    % Loop through boundaries
    
    
    
    
    

end