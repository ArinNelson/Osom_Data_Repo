%=========================================================================%
% rivers_to_roms_grid.m
% Determine river locations on ROMS grid
% 
% Arin Nelson
% on 07/29/2021
% 
% Last edited 07/29/2021
%=========================================================================%
clc; close all; addpath('D:\OSOM_Data_Repo\git\MATLAB\Utilities');

% Switches
Switch = zeros(9,1);
Switch(1) = 0;  % Load data
Switch(2) = 1;  % Lookit


% Options
rvrgg_dir = 'D:/OSOM_Data_Repo/USGS/RiverGauges/';
grid_file = 'D:/ROMS/Resources/ngbay_grd.nc';

%=========================================================================%
if(Switch(1))
    
	% List of stations
    rvrgg_list = ls([rvrgg_dir '/0*']);
    n_rvrgg    = size(rvrgg_list,1);
    
    % Load info from each station
    rvrgg = struct;
    for i=1:n_rvrgg
        tmp = load([rvrgg_dir '/' rvrgg_list(i,:) '/' rvrgg_list(i,:) '_info.mat']);
        rvrgg(i).id   = rvrgg_list(i,:);
        rvrgg(i).name = tmp.station_name;
        rvrgg(i).lon  = tmp.longitude;
        rvrgg(i).lat  = tmp.latitude;
        rvrgg(i).alt  = tmp.altitude;
        rvrgg(i).area = tmp.drainage_area;
    end
    clear i tmp;
    
    % Load grid data
    grid = grid_get(grid_file);
    
    % Index meshgrid
    ii = 1:size(grid.lon_rho,1);    
    jj = 1:size(grid.lon_rho,2);    
    [jj,ii] = meshgrid(jj,ii);
    iu = ii(1:end-1,:);     iv = ii(:,1:end-1);
    ju = jj(1:end-1,:);     jv = jj(:,1:end-1);
    
    % Find valid points
    m_rho = find(grid.mask_rho==1 & grid.lon_rho>-71.9 & grid.lon_rho<-71 & grid.lat_rho>41.25);
    m_u   = find(grid.mask_u  ==1 & grid.lon_u  >-71.9 & grid.lon_u  <-71 & grid.lat_u  >41.25);
    m_v   = find(grid.mask_v  ==1 & grid.lon_v  >-71.9 & grid.lon_v  <-71 & grid.lat_v  >41.25);

    % Keep valid points
    x_rho = grid.lon_rho(m_rho);    x_u = grid.lon_u(m_u);      x_v = grid.lon_v(m_v);
    y_rho = grid.lat_rho(m_rho);    y_u = grid.lat_u(m_u);      y_v = grid.lat_v(m_v);
    i_rho = ii(m_rho)+i.*jj(m_rho); i_u = iu(m_u)+i.*ju(m_u);   i_v = iv(m_v)+i.*jv(m_v);

    % Interpolate river lon/lat to grid i/j
    ntrplnt = scatteredInterpolant(grid.lon_rho(:),grid.lat_rho(:),ii(:)+i.*jj(:),'linear');
    i_rvr = zeros(n_rvrgg,1);
    j_rvr = zeros(n_rvrgg,1);
    for i=1:n_rvrgg
        tmp = ntrplnt(rvrgg(i).lon,rvrgg(i).lat);
        i_rvr(i) = real(tmp);
        j_rvr(i) = imag(tmp);
    end
    clear i tmp;
    
    % Clean-up
    clear ii jj iu ju iv jv m_rho m_u m_v grid
    
end
%=========================================================================%
if(Switch(2))
    
    % Plot
    figure('units','normalized','position',[0 0 1 1]);
    axes('units','normalized','position',[0 0 1 0.975]);
    hold on;
%         scatter(x_rho,y_rho,12,i_rho,'sk','filled');
%         scatter(x_u,y_u,12,i_u,'ob','filled');
%         scatter(x_v,y_v,12,i_v,'or','filled');
%         for i=1:numel(rvrgg)
%             plot(rvrgg(i).lon,rvrgg(i).lat,'ok');
%         end
        scatter( real(i_rho), imag(i_rho), 12, 'sk' , 'filled' );
        scatter( real(i_u)+0.5,   imag(i_u),   12, 'ob' , 'filled' );
        scatter( real(i_v), imag(i_v)+0.5, 12, 'or', 'filled' );
        scatter( i_rvr,     j_rvr, 12, '+k' );
        for i=1:n_rvrgg
            text(i_rvr(i),j_rvr(i)+0.5,rvrgg(i).id);
        end
    hold off; box on; axis equal; axis tight;
    
end