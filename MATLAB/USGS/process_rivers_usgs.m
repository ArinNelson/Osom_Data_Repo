%=========================================================================%
% process_rivers_usgs.m
% Gather data from multiple rivers
% 
% Arin Nelson
% on 07/29/2021
% 
% Last edited 07/29/2021
%=========================================================================%

% Options
rvrgg_year = [1987 2020];     % Collect data through these years
rvrgg_dir  = 'D:/OSOM_Data_Repo/USGS/RiverGauges/';

% Variable info
cfs_to_cms = 0.028316846592;    % convert cubic feet per second to cubic meters per second
ft_to_m    = 1/3.281;           % convert feet to meters
var_info   = {'00060',  'q',   cfs_to_cms; ...     % m3/s
              '00065',  'h',   ft_to_m;    ...     % m
              '00010',  'T',   1;          ...     % C
              '00095',  'C',   1;          ...     % ?
              '00300',  'DO2', 1;          ...     % mg/l
              '00400',  'pH',  1;          ...     % ?
             };

%=========================================================================%

	% Stations available
    %rvrgg_list = ls([rvrgg_dir '/0*']);
    rvrgg_list = ['01108000';...
                  '01109060';...
                  '01109220';...
                  '01109403';...
                  '01112500';...
                  '01113895';...
                  '01114000';...
                  '01114500';...
                  '01116500';...
                  '01117000';...
                  '01118500';...
                 ];
    n_rvrgg    = size(rvrgg_list,1);
    stats      = cell(n_rvrgg,1);
    
    % Create master files
    for i=1:n_rvrgg
        stats{i} = process_river_usgs(rvrgg_list(i,:),rvrgg_year,rvrgg_dir);
    end
    
    % Combo stat's
    yy = rvrgg_year(1) : 1 : rvrgg_year(end);
    ny = numel(yy);
    qmean = NaN(ny,n_rvrgg);
    qvar  = NaN(ny,n_rvrgg);
    for i=1:ny
    for j=1:n_rvrgg
        qmean(i,j) = stats{j}(i,1,2);
        qvar(i,j)  = stats{j}(i,1,3);
    end
    end
    clear i j;
    %errorbar(yy,qmean,sqrt(qvar));
    plot(yy,log(qmean),'-'); legend(rvrgg_list);

%=========================================================================%