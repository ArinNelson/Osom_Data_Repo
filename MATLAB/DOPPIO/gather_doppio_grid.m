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
max_wait        = 500;	% max time to wait (secs) for web response (uses Java & parallelism, set to <=0 to disable)
dopp_gridfile   = 'F:/OSOM_Data_Repo/DOPPIO/doppio_grid.nc';      % save file name
osom_gridfile   = 'C:/Library/ROMS_Stuff/Resources/ngbay_grd.nc';
osom_vertfile   = 'C:/Library/ROMS_Stuff/Resources/forcefiles_riroms/riroms_ini_20060101.nc';   % File with vertical grid info

% Constants
dopp_url = 'https://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best';

% Gather grids
grid_osom   = grid_get(osom_gridfile);
%grid_dopp   = grid_get(dopp_url );
load('grid_dopp.mat');

% Determine DOPPIO grid within OSOM domain
osom_xpoly	= [grid_osom.lon_rho(1,:)'; grid_osom.lon_rho(2:end-1,end); grid_osom.lon_rho(end,end:-1:1)'; grid_osom.lon_rho(end:-1:1,1)];
osom_ypoly	= [grid_osom.lat_rho(1,:)'; grid_osom.lat_rho(2:end-1,end); grid_osom.lat_rho(end,end:-1:1)'; grid_osom.lat_rho(end:-1:1,1)]; 
[ii,jj]     = find( inpolygon(grid_dopp.lon_rho,grid_dopp.lat_rho,osom_xpoly,osom_ypoly)==1 );
ii = min(ii)-1 : max(ii)+1;     i0 = ii(1);     ni = numel(ii);
jj = min(jj)-1 : max(jj)+1;     j0 = jj(1);     nj = numel(jj);

% Truncate DOPPIO grid
grid_dopp   = grid_truncate(grid_dopp,ii,jj,[]);

% Save grid
grid_save(dopp_gridfile,grid_dopp);
