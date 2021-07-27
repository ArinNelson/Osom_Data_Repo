%=========================================================================%
% gather_doppio_best.m
% Gather DOPPIO 'best' data and sort by variable and year/month
% by Arin Nelson
% on 07/26/2021
% 
% last edited by Arin Nelson on 07/26/2021
%=========================================================================%
clc; clear mex; addpath('../Utilities');

% Options
date_start = [2017,12];     % Start year, month to gather data for
date_end   = [2021,01];     % End   year, month to gather data for
max_wait   = 300;           % Max time to wait (secs) for web response
%save_dir   = 'D:/OSOM_Data_Repo/DOPPIO/best/';
save_dir   = 'F:/OSOM_Data_Repo/DOPPIO/best/';
var_to_get = {'zeta','ubar','vbar','u','v','temp','salt'};
%osom_grid  = 'D:/ROMS/Resources/ngbay_grd.nc';	% OSOM grid file

% Constants
dopp_url   = 'https://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best';
n_var      = numel(var_to_get);

% Variable info     (name,  # dims,     grid type)
var_info = {'zeta',2,1; ...
            'ubar',2,2; ...
            'vbar',2,3; ...
            'u',   3,2; ...
            'v',   3,3; ...
            'temp',3,1; ...
            'salt',3,1; ...
           };
       
% NOTES
% time is in units of hours since 2011/01/01 00:00:00

%=========================================================================%

% Ensure data is accessible
dopp_info = ncinfo_web(max_wait,dopp_url);
nz = dopp_info.Dimensions( strcmp({dopp_info.Dimensions.Name},'s_rho')==1 ).Length;

% Ensure variables in var_to_get are valid
i_var = zeros(n_var,1);
for iv=1:n_var
    ii = find( strcmp(var_info(:,1),var_to_get{iv})==1 );
    if(~isempty(ii))
        i_var(iv) = ii;
    else
        error(['Invalid or unimplemented variable: ' var_to_get{iv}]);
    end
end
clear iv ii;

% Ensure data directories exist
if(exist(save_dir,'dir')~=7);   mkdir(save_dir);    end
for iv=1:n_var
    if(exist([save_dir var_to_get{iv}],'dir')~=7);  mkdir([save_dir var_to_get{iv}]);   end
end
clear iv;

% Determine DOPPIO grid within OSOM domain
if(exist('dopp_gatherinfo.mat','file')~=2)
    grid_osom = grid_get(osom_grid);
    grid_dopp = grid_get(dopp_url );
    osom_xpoly	= [grid_osom(1).lon(1,:)'; grid_osom(1).lon(2:end-1,end); grid_osom(1).lon(end,end:-1:1)'; grid_osom(1).lon(end:-1:1,1)];
    osom_ypoly	= [grid_osom(1).lat(1,:)'; grid_osom(1).lat(2:end-1,end); grid_osom(1).lat(end,end:-1:1)'; grid_osom(1).lat(end:-1:1,1)]; 
    [ix,iy]     = find( inpolygon(grid_dopp(1).lon,grid_dopp(1).lat,osom_xpoly,osom_ypoly)==1 );
    ii = min(ix)-1 : max(ix)+1;     i0 = ii(1);     ni = numel(ii);
    jj = min(iy)-1 : max(iy)+1;     j0 = jj(1);     nj = numel(jj);
    grid_dopp = grid_truncate(grid_dopp,ii,jj);
    clear grid_osom osom_xpoly osom_ypoly ix iy ii jj;
    save('dopp_gatherinfo.mat','grid_dopp','i0','j0','ni','nj');
else
    load('dopp_gatherinfo.mat');
end

% Load time variable
time = ncread_web(max_wait,dopp_url,'time');
time = (time./24) + datenum(2017,11,01,0,0,0);

% Clean-up
clear dopp_info grid_osom tmp osom_xpoly osom_ypoly ix iy ii jj;

%=========================================================================%

% Loop through years and months
year_on  = date_start(1);
month_on = date_start(2);
while( datenum(year_on,month_on,1) <= datenum(date_end(1),date_end(2),1) )
    
    % Time indices for this month
    it  = find(year(time)==year_on & month(time)==month_on);
    nt  = numel(it);
    it0 = it(1);
    
    % Loop through variables
    for iv=1:n_var
    
        % This year/month save file
        save_file = [save_dir var_to_get{iv} '/DOPPIO_best_' var_to_get{iv} '_' num2str(year_on) '_' sprintf('%0.2d',month_on) '.nc'];
        if(exist(save_file,'file')~=2)
        disp(['Gathering ' var_to_get{iv} ' for ' num2str(year_on) '/' sprintf('%0.2d',month_on) '...']);
       
            % Start array
            i_start = ones( var_info{i_var(iv),2}+1, 1 );   
            i_start(1:2) = [i0, j0];
            i_start(end) = it0;
            
            % Count array
            i_count = inf( var_info{i_var(iv),2}+1, 1 );
            switch var_info{i_var(iv),3}
                case 1; i_count(1:2) = [ni,   nj  ];
                case 2; i_count(1:2) = [ni-1, nj  ];
                case 3; i_count(1:2) = [ni,   nj-1];
            end
            i_count(end) = nt;
            
            % Gather variable for this year-month
            while(exist('dat','var')~=1)
                try dat = ncread_web(max_wait,dopp_url,var_to_get{iv},i_start,i_count);
                catch err
                end
            end
            
            % Save to save_file
            nc_gen_dopp_data(save_file,ni,nj,nz,var_to_get(iv));    % NOTE PARENTHESIS USED, NOT {}!
            ncwrite(save_file,'time',        time(it));
            ncwrite(save_file,var_to_get{iv},dat     );
            
            % Clean-up
            clear i_start i_count dat;
            
        end
        clear save_file;
    
    end
    clear iv;
    
    % Clean-up
    clear it nt it0;
    
	% Onto next year/month
	month_on = month_on + 1;
    if(month_on == 13)
        year_on = year_on + 1;  
        month_on = 1;   
    end
    
end
clear year_on month_on;

%=========================================================================%