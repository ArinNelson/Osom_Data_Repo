function [ndx_i, ndx_j, ndx_w] = grid_interpolant(x,y,xq,yq,varargin)

    % If either grid is not meshgridded, meshgrid it
    if(numel(size(x))==2 && any(size(x)==1))
        [y,x] = meshgrid(y,x);
    end
    if(numel(size(xq))==2 && any(size(xq)==1))
        [yq,xq] = meshgrid(yq,xq);
    end
    
    % Use masks if employed
    if(nargin>4)
        m = varargin{1};
        mq = varargin{2};
    else
        m  = ones(size(x ));
        mq = ones(size(xq));
    end

    % For each point in xq, yq, find 2 surrounding points in x, y
    ndx_i = zeros(size(xq,1),size(xq,2),2);
    ndx_j = zeros(size(xq,1),size(xq,2),2);
    ndx_w = zeros(size(xq,1),size(xq,2),2,2);
    for i=1:size(xq,1)
    for j=1:size(xq,2)
    if(mq(i,j)==1)    
        
        % Find nearest point
        dd = sqrt( (x-xq(i,j)).^2 + (y-yq(i,j)).^2 );
        [ii,jj] = find(dd == min(dd(:)));
        
        % Find 3 other points in opposite directions 
        if(~isempty(ii))
            
            % Gather the 4 surrounding points
            if(     x(ii,jj) <xq(i,j) && y(ii,jj)>=yq(i,j));	ii = [ii   ii+1; ii     ii+1 ];      jj = [jj   jj  ; jj-1 jj-1];
            elseif( x(ii,jj)>=xq(i,j) && y(ii,jj)>=yq(i,j));    ii = [ii-1 ii  ; ii-1   ii    ];     jj = [jj   jj  ; jj-1 jj-1];
            elseif( x(ii,jj)>=xq(i,j) && y(ii,jj) <yq(i,j));    ii = [ii-1 ii  ; ii-1   ii    ];     jj = [jj+1 jj+1; jj   jj  ];
            elseif( x(ii,jj) <xq(i,j) && y(ii,jj) <yq(i,j));    ii = [ii   ii+1; ii     ii+1  ];     jj = [jj+1 jj+1; jj   jj  ];
            end
        
            % Checks
            if(any(ii(:)==0));          ii = ii + 1;    end
            if(any(ii(:)>size(x,1)));   ii = ii-1;      end
            if(any(jj(:)==0));          jj = jj + 1;    end
            if(any(jj(:)>size(x,2)));   jj = jj-1;      end
            
            % (checks for being on edge later...)
            
            % Compute relative weights
            ndx_i(i,j,:) = [nanmin(ii(:)) nanmax(ii(:))];
            ndx_j(i,j,:) = [nanmin(jj(:)) nanmax(jj(:))];
            for a=1:2
            for b=1:2
            if(m(ii(a,b),jj(a,b))==1)    
                ndx_w(i,j,a,b) = 1 / sqrt( (x(ii(a,b),jj(a,b))-xq(i,j))^2 + (y(ii(a,b),jj(a,b))-yq(i,j))^2 );
            end
            end
            end
            clear a b;
            
            % Ensure sum of weights is 1
            tmp = squeeze(ndx_w(i,j,:,:));
            tmp = nansum(tmp(:));
            ndx_w(i,j,:,:) = ndx_w(i,j,:,:) ./ tmp;
            clear tmp;
            
        end
        clear dd ii jj;
     
    end
    end
    end
    clear i j;

end