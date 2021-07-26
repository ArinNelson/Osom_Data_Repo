%=========================================================================%
% gather_doppio_best.m
% Gather DOPPIO forecast 'best' data, to use for model spin-up
% by Arin Nelson
% on 07/22/2021
%
% last edited by Arin Nelson on 07/22/2021
%=========================================================================%
clc; addpath('../Utilities');   % clear mex; % (useful when debugging nc_gen_* codes)

% Options
max_wait   = 120;           % Max time to wait (secs) for web response (uses Java & parallelism, set to 0 to disable)
%data_dir   = '/gpfs/data/epscor/anelson5/OSOM_Data_Repo/DOPPIO/best/';
data_dir   = 'F:/OSOM_Data_Repo/DOPPIO/best/';          % My desktop
grid_file  = 'F:/OSOM_Data_Repo/DOPPIO/doppio_grid.nc';	% Grid file
srf_file   = 'OSOM_frc_DOPPIO_best.nc';
bry_file   = 'OSOM_bry_DOPPIO_best.nc';
ini_file   = 'OSOM_ini_DOPPIO_best.nc'; 

% Constants
date_min   = [2017,11,01];
%date_max   = [2021,07,01];
dopp_url   = 'https://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_fmrc.ncd';   % Contains both 'best' and forecasts in 2 time dimensions

% % Info on possible variables to gather
% %         variable, # dims, horiz. grid type (1=rho, 2=u, 3=v), vert. grid type (1=rho, 2=w
% var_info = {'zeta',        2, 1, 0; ...
%             'ubar',        2, 2, 0; ...
%             'vbar',        2, 3, 0; ...
%             'u',           3, 2, 1; ...
%             'v',           3, 3, 1; ...
%             'w',           3, 1, 2; ...
%             'temp',        3, 1, 1; ...
%             'salt',        3, 1, 1; ...
%             'AKv',         3, 1, 2; ...
%             'AKt',         3, 1, 2; ...
%             'AKs',         3, 1, 2; ...
%             'shflux',      2, 1, 0; ...
%             'ssflux',      2, 1, 0; ...
%             'swrad_daily', 2, 1, 0; ...
%             'sustr',       2, 2, 0; ...
%             'svstr',       2, 3, 0; ...
%            };
% 
% % Save file info
% file_info = {'srf',{'shflux','ssflux','swrad_daily','sustr','svstr'};   ...
%              'bry',{'zeta','ubar','vbar','u','v','temp','salt'};        ...
%              'ini',{'zeta','ubar','vbar','u','v','temp','salt'};        ...
%             };
       
%=========================================================================%

% if grid file does not exist, ask user to do it
% if it does exist, load the grid variables
%if(exist(grid_file,'file')~=2)
%   error('Specified grid file does not exist.  Either fix the grid_file variable or run gather_nam_grid.m.');
%else
%   (load grid vars)
%end
load('bkup_grid.mat');

% % Save directories
% for iv=1:numel(var_to_get)
% 	var_dir = [data_dir var_to_get{iv}];
% 	if(exist(var_dir,'dir')~=7)
%         mkdir(var_dir);
% 	end
% end
% clear iv var_dir;

%-------------------------------------------------------------------------%

% Load the time variable
time = ncread(base_url,'time');
time = (time./24) + datenum(2017,11,01);

% Initialize files




