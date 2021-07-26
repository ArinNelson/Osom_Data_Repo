function [ndx_i, ndx_j, ndx_w] = grid_intersect(grid1, grid2)

    % Init structures
    ndx_i = cell(3,1);
    ndx_j = cell(3,1);
    ndx_w = cell(3,1);
  
    % Loop through grids
    for ig=1:3

        % Init variables
        ndx_i{ig} = zeros(size(grid1(ig).lon,1), size(grid1(ig).lon,2), 2);
        ndx_j{ig} = zeros(size(grid1(ig).lon,1), size(grid1(ig).lon,2), 2);
        ndx_w{ig} = zeros(size(grid1(ig).lon,1), size(grid1(ig).lon,2), 2, 2);
        
        % Loop through
        for i=1:size(grid2(ig).lon,1)-1
        for j=1:size(grid2(ig).lon,2)-1

            % This polygon from grid2 pts
            x_poly = [grid2(ig).lon(i,j), grid2(ig).lon(i+1,j), grid2(ig).lon(i+1,j+1), grid2(ig).lon(i,j+1)];
            y_poly = [grid2(ig).lat(i,j), grid2(ig).lat(i+1,j), grid2(ig).lat(i+1,j+1), grid2(ig).lat(i,j+1)];
  
            % Find grid1 points within this grid2 polygon
            in_poly = inpolygon(grid1(ig).lon,grid1(ig).lat,x_poly,y_poly);
            
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
        for i=1:size(grid1(ig).lon,1)
        for j=1:size(grid1(ig).lon,2)
        if(grid1(ig).mask(i,j)==1)
            
            % The 4 surrounding points
            for a=1:2
            for b=1:2
                aa = ndx_i{ig}(i,j,a);
                bb = ndx_j{ig}(i,j,b);
                if(grid2(ig).mask(aa,bb)==1)
                    ndx_w{ig}(i,j,a,b) = circledist( grid1(ig).lon(i,j), grid1(ig).lat(i,j), grid2(ig).lon(aa,bb), grid2(ig).lat(aa,bb) );
                end                              
            end
            end
            clear a b aa bb;
            
            % If all surrounding points were land, find nearest water point
            test = ndx_w{ig}(i,j,:,:);
            if(all(test(:)==0))
                kk      = find(grid2(ig).mask == 1);
                [ii,jj] = find(grid2(ig).mask == 1); 
                dd = zeros(numel(kk));
                for k=1:numel(kk);  dd(k) = circledist( grid1(ig).lon(i,j), grid1(ig).lat(i,j), grid2(ig).lon(kk(k)), grid2(ig).lat(kk(k)) );   end
                mm = find( dd==min(dd), 1, 'first' );
                ndx_i{ig}(i,j,:)   = ones(2,1).*ii(mm);
                ndx_j{ig}(i,j,:)   = ones(2,1).*jj(mm);
                ndx_w{ig}(i,j,:,:) = ones(2,2);
                clear kk ii jj dd mm;
            end
            clear test;
          
        end
        end
        end
        clear i j;

    end
    clear ig;
      
end