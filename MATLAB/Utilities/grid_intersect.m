function [ndx_i, ndx_j, ndx_w] = grid_intersect(grid1, grid2)

    % Init structures
    ndx_i = cell(3,1);
    ndx_j = cell(3,1);
    ndx_w = cell(3,1);
  
    % Loop through grids
    str_grid = {'rho','u','v'};
    for ig=1:numel(str_grid)

        % Grid vars
        eval(['lon1  = grid1.lon_'  str_grid{ig} ';']);
        eval(['lat1  = grid1.lat_'  str_grid{ig} ';']);
        eval(['mask1 = grid1.mask_' str_grid{ig} ';']);
        ni1 = size(lon1,1);  nj1 = size(lon1,2);
        
        % Init variables
        ndx_i{ig} = zeros(ni1,nj1,2);
        ndx_j{ig} = zeros(ni1,nj1,2);
        ndx_w{ig} = zeros(ni1,nj1,2,2);
        
        % Other grid vars
        eval(['lon2  = grid2.lon_'  str_grid{ig} ';']);
        eval(['lat2  = grid2.lat_'  str_grid{ig} ';']);
        eval(['mask2 = grid2.mask_' str_grid{ig} ';']);
        ni2 = size(lon2,1); nj2 = size(lon2,2);
        
        % Loop through
        for i=1:ni2-1
        for j=1:nj2-1

            % This polygon from grid2 pts
            x_poly = [lon2(i,j), lon2(i+1,j), lon2(i+1,j+1), lon2(i,j+1)];
            y_poly = [lat2(i,j), lat2(i+1,j), lat2(i+1,j+1), lat2(i,j+1)];
  
            % Find grid1 points within this grid2 polygon
            in_poly = inpolygon(lon1,lat1,x_poly,y_poly);
            
            % Save points that are not masked as land
            if(any(in_poly(:)==1))
                [i_in,j_in] = find(in_poly==1);
                for k=1:numel(i_in)
                    ndx_i{ig}(i_in(k),j_in(k),:) = [i i+1];
                    ndx_j{ig}(i_in(k),j_in(k),:) = [j j+1];
                end
            end
            
            % Clean-up
            clear xpoly ypoly inpoly;

        end
        end
        clear i j;
        
        % Compute interpolation weights
        if(nargout>2)
        for i=1:ni1
        for j=1:nj1
        if(mask1(i,j)==1)
            
            % The 4 surrounding points
            for a=1:2
            for b=1:2
                aa = ndx_i{ig}(i,j,a);
                bb = ndx_j{ig}(i,j,b);
                if(mask2(aa,bb)==1)
                    ndx_w{ig}(i,j,a,b) = 1/circledist( lon1(i,j), lat1(i,j), lon2(aa,bb), lat2(aa,bb) );
                end                              
            end
            end
            clear a b aa bb;
            
            % Set sum of weights to be 1
            tmp = squeeze( ndx_w{ig}(i,j,:,:) );
            ndx_w{ig}(i,j,:,:) = ndx_w{ig}(i,j,:,:) ./ sum(tmp(:));
            clear tmp;
            
            % If all surrounding points were land, find nearest water point
            test = ndx_w{ig}(i,j,:,:);
            if(all(test(:)==0))
                kk      = find(mask2 == 1);
                [ii,jj] = find(mask2 == 1); 
                dd = zeros(numel(kk));
                for k=1:numel(kk)
                    
                    %dd(k) = circledist( lon1(i,j), lat1(i,j), lon2(kk(k)), lat2(kk(k)) );   
                    dd(k) = (
                
                end
                mm = find( dd==min(dd), 1, 'first' );
                ndx_i{ig}(i,j,:)   = ones(2,1).*ii(mm);
                ndx_j{ig}(i,j,:)   = ones(2,1).*jj(mm);
                ndx_w{ig}(i,j,:,:) = ones(2,2)./4;
                clear kk ii jj dd mm;
            end
            clear test;
          
        end
        end
        end
        clear i j;
        end

    end
    clear ig;
    
end