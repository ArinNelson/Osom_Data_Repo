%=========================================================================%
% gather_doppio_grid.m
% Gather DOPPIO grid 
% by Arin Nelson
% on 07/22/2021
%
% last edited by Arin Nelson on 07/22/2021
%=========================================================================%
clc; addpath('../Utilities');   % clear mex; % (useful when debugging nc_gen_* codes)

% Options
max_wait        = 120;	% max time to wait (secs) for web response (uses Java & parallelism, set to <=0 to disable)
grid_savedir    = 'F:/OSOM_Data_Repo/DOPPIO/';      % location to save grid file
dopp_gridfile   = [grid_savedir 'doppio_grid.nc'];      % grid file name
osom_gridfile   = 'C:\Library\ROMS_Stuff\Resources\ngbay_grd.nc';
osom_vertfile   = 'C:\Library\ROMS_Stuff\Resources\forcefiles_riroms\riroms_ini_20060101.nc';   % File with vertical grid info

% Constants
grid_url = 'https://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/runs/History_RUN_2017-11-01T00:00:00Z';

%-------------------------------------------------------------------------%

% Full grid file directory
grid_file = [grid_savedir dopp_gridfile];

% if the grid file already exists, ask user if they want to ovewrite it
check = false;
if(exist(grid_file,'file')==2)
	usr_input = input('A grid file with this name at this directory already exists.  Overwrite it? 1=yes, 0=no');
    if(usr_input==1)
        delete(grid_file);
    else
        check = true;
    end
end
if(check);  error('User aborted.'); end

% Read in RHO-grid variables from url
dopp_lon_rho = ncread_web(max_wait,grid_url,'lon_rho');
dopp_lat_rho = ncread_web(max_wait,grid_url,'lat_rho');

% OSOM grid variables
osom_lon_rho   = ncread(osom_gridfile,'lon_rho' );
osom_lat_rho   = ncread(osom_gridfile,'lat_rho' );
osom_mask_rho  = ncread(osom_gridfile,'mask_rho');
osom_angle_rho = ncread(osom_gridfile,'mask_rho');
osom_lon_u     = ncread(osom_gridfile,'lon_u'   ); 
osom_lat_u     = ncread(osom_gridfile,'lat_u'   );
osom_mask_u    = ncread(osom_gridfile,'mask_u'  );
osom_lon_v     = ncread(osom_gridfile,'lon_v'   );  
osom_lat_v     = ncread(osom_gridfile,'lat_v'   );
osom_mask_v    = ncread(osom_gridfile,'mask_v'  );
osom_s_rho     = ncread(osom_vertfile,'s_rho'   );

% Find indices of grid overlapping given domain
osom_xpoly = [osom_lon_rho(1,:)'; osom_lon_rho(2:end-1,end); osom_lon_rho(end,end:-1:1)'; osom_lon_rho(end:-1:1,1)];
osom_ypoly = [osom_lat_rho(1,:)'; osom_lat_rho(2:end-1,end); osom_lat_rho(end,end:-1:1)'; osom_lat_rho(end:-1:1,1)]; 
[ix,iy] = find( inpolygon(dopp_lon_rho,dopp_lat_rho,osom_xpoly,osom_ypoly)==1 );
ii = (min(ix)-1) : (max(ix)+1);   i0 = ii(1);    ni = numel(ii);
jj = (min(iy)-1) : (max(iy)+1);   j0 = jj(1);    nj = numel(jj);
clear osom_xpoly osom_ypoly ix iy;

% Gather remaining grid variables
dopp_lon_rho   = dopp_lon_rho(ii,jj);
dopp_lat_rho   = dopp_lat_rho(ii,jj);
dopp_mask_rho  = ncread_web(max_wait,grid_url,'mask_rho',[i0 j0],[ni   nj  ]);
dopp_pm        = ncread_web(max_wait,grid_url,'pm',      [i0 j0],[ni   nj  ]);
dopp_pn        = ncread_web(max_wait,grid_url,'pn',      [i0 j0],[ni   nj  ]);
dopp_angle_rho = ncread_web(max_wait,grid_url,'angle',   [i0 j0],[ni   nj  ]);
dopp_lon_u     = ncread_web(max_wait,grid_url,'lon_u',   [i0 j0],[ni-1 nj  ]);
dopp_lat_u     = ncread_web(max_wait,grid_url,'lat_u',   [i0 j0],[ni-1 nj  ]);
dopp_mask_u    = ncread_web(max_wait,grid_url,'mask_u',  [i0 j0],[ni-1 nj  ]);
dopp_lon_v     = ncread_web(max_wait,grid_url,'lon_v',   [i0 j0],[ni   nj-1]);
dopp_lat_v     = ncread_web(max_wait,grid_url,'lat_v',   [i0 j0],[ni   nj-1]);
dopp_mask_v    = ncread_web(max_wait,grid_url,'mask_v',  [i0 j0],[ni   nj-1]);
dopp_s_rho     = ncread_web(max_wait,grid_url,'s_rho');

% Determine angles on u and v grids
osom_angle_u = ( osom_angle_rho(1:end-1,:) + osom_angle_rho(2:end,:) )./2;
osom_angle_v = ( osom_angle_rho(:,1:end-1) + osom_angle_rho(:,2:end) )./2;
dopp_angle_u = ( dopp_angle_rho(1:end-1,:) + dopp_angle_rho(2:end,:) )./2;
dopp_angle_v = ( dopp_angle_rho(:,1:end-1) + dopp_angle_rho(:,2:end) )./2;

% Create weights for linear interpolation from DOPPIO grid to ROMS-OSOM grid
ndxi_rho = NaN(size(osom_lon_rho,1),size(osom_lon_rho,2),2);
ndxj_rho = NaN(size(osom_lon_rho,1),size(osom_lon_rho,2),2);
ndxi_u   = NaN(size(osom_lon_u,1),  size(osom_lon_u,2),  2);
ndxj_u   = NaN(size(osom_lon_u,1),  size(osom_lon_u,2),  2);
ndxi_v   = NaN(size(osom_lon_v,1),  size(osom_lon_v,2),  2);
ndxj_v   = NaN(size(osom_lon_v,1),  size(osom_lon_v,2),  2);
for i=1:size(dopp_lon_rho,1)-1
for j=1:size(dopp_lon_rho,2)-1
    
	% Polygon/box from these 4 DOPPIO RHO-grid points
	dopp_xpoly = [dopp_lon_rho(i,j), dopp_lon_rho(i+1,j), dopp_lon_rho(i+1,j+1), dopp_lon_rho(i,j+1)];
	dopp_ypoly = [dopp_lat_rho(i,j), dopp_lat_rho(i+1,j), dopp_lat_rho(i+1,j+1), dopp_lat_rho(i,j+1)];
  
	% Find OSOM grid points within this polygon
	osom_inpoly = inpolygon(osom_lon_rho,osom_lat_rho,dopp_xpoly,dopp_ypoly);
  
	% If any, save their indices
	if(any(osom_inpoly(:)==1))
        [i_in,j_in] = find(osom_inpoly==1);
        for k=1:numel(i_in)
            ndxi_rho(i_in(k),j_in(k),:) = [i i+1];
            ndxj_rho(i_in(k),j_in(k),:) = [j j+1];
        end
	end
    clear dopp_xpoly dopp_ypoly osom_inpoly i_in j_in k;
    
    % Same for u
    if(i<size(dopp_lon_rho,1)-1)
        dopp_xpoly = [dopp_lon_u(i,j), dopp_lon_u(i+1,j), dopp_lon_u(i+1,j+1), dopp_lon_u(i,j+1)];
        dopp_ypoly = [dopp_lat_u(i,j), dopp_lat_u(i+1,j), dopp_lat_u(i+1,j+1), dopp_lat_u(i,j+1)];
        osom_inpoly = inpolygon(osom_lon_u, osom_lat_u, dopp_xpoly, dopp_ypoly);
        if(any(osom_inpoly(:)==1))
            [i_in,j_in] = find(osom_inpoly==1);
            for k=1:numel(i_in)
                ndxi_u(i_in(k),j_in(k),:) = [i i+1];
                ndxj_u(i_in(k),j_in(k),:) = [j j+1];
            end
        end
        clear dopp_xpoly dopp_ypoly osom_inpoly i_in j_in k;
    end

    % Same for v
    if(j<size(dopp_lon_rho,2)-1)
        dopp_xpoly = [dopp_lon_v(i,j), dopp_lon_v(i+1,j), dopp_lon_v(i+1,j+1), dopp_lon_v(i,j+1)];
        dopp_ypoly = [dopp_lat_v(i,j), dopp_lat_v(i+1,j), dopp_lat_v(i+1,j+1), dopp_lat_v(i,j+1)];
        osom_inpoly = inpolygon(osom_lon_v, osom_lat_v, dopp_xpoly, dopp_ypoly);
        if(any(osom_inpoly(:)==1))
            [i_in,j_in] = find(osom_inpoly==1);
            for k=1:numel(i_in)
                ndxi_v(i_in(k),j_in(k),:) = [i i+1];
                ndxj_v(i_in(k),j_in(k),:) = [j j+1];
            end
        end
        clear dopp_xpoly dopp_ypoly osom_inpoly i_in j_in k;
    end
  
end
end
clear i j;

% Compute the linear interpolation weights using circle distance
ntrp_wght_rho = NaN(size(osom_lon_rho,1),size(osom_lon_rho,2),2,2);
ntrp_wght_u   = NaN(size(osom_lon_u,  1),size(osom_lon_u,  2),2,2);
ntrp_wght_v   = NaN(size(osom_lon_v,  1),size(osom_lon_v,  2),2,2);
for i=1:size(osom_lon_rho,1)
for j=1:size(osom_lon_rho,2)
    
	% Rho interpolation weights
    for a=1:2
    for b=1:2
        ntrp_wght_rho(i,j,a,b) = circledist(osom_lon_rho(i,j), osom_lat_rho(i,j), ...
                                            dopp_lon_rho(ndxi_rho(i,j,a),ndxj_rho(i,j,b)), ...
                                            dopp_lat_rho(ndxi_rho(i,j,a),ndxj_rho(i,j,b)));
    end
    end
    
    % Same for u,v
    if(i<=size(osom_lon_u,1))
    for a=1:2
    for b=1:2
        ntrp_wght_u(i,j,a,b) = circledist(osom_lon_u(i,j), osom_lat_u(i,j), ...
                                          dopp_lon_u(ndxi_u(i,j,a),ndxj_u(i,j,b)), ...
                                          dopp_lat_u(ndxi_u(i,j,a),ndxj_u(i,j,b)));
    end
    end
    end
    
    % Same for v
    if(j<=size(osom_lon_v,2))
    for a=1:2
    for b=1:2
        ntrp_wght_v(i,j,a,b) = circledist(osom_lon_v(i,j), osom_lat_v(i,j), ...
                                          dopp_lon_v(ndxi_v(i,j,a),ndxj_v(i,j,b)), ...
                                          dopp_lat_v(ndxi_v(i,j,a),ndxj_v(i,j,b)));
    end
    end
    end
    
end
end
clear i j a b;

%-------------------------------------------------------------------------%

% For now, save as .mat file
save('bkup_grid.mat','osom_*','dopp_*','ntrp_wght_*','ndxi_*','ndxj_*','i0','j0','ni','nj');

% Generate the grid file
%nc_gen_doppio_grid(grid_file,


% Save the grid file