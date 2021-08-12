%=========================================================================%
% osom_gen_secs.m
% Generate OSOM domain sections
% by Arin Nelson
% on 08/10/2021
% 
% Last updated by Arin Nelson on 08/10/2021
%=========================================================================%
addpath('../Utilities');

% Options
baysecs_shapefile = 'D:\OSOM_Data_Repo\GIS\NBEP_BaySecs\BAYSECTIONS_NBEP2017.shp';
baysegs_shapefile = 'D:\OSOM_Data_Repo\GIS\NBEP_BaySegs\BAYSEGMENTS_NBEP2017.shp';
osom_gridfile     = 'D:\ROMS\Resources\ngbay_grd.nc';

% Secions and segments wanted
sec_to_get = {'Providence R.',      6,          'sec';
              'Taunton R.',         8,          'sec';
              'Upper Narr. Bay',	9,          'sec';
              'Mount Hope Bay',     3,          'sec';
              'Greenwich Bay',      2,          'sec';
              'Sakonnet River',     7,          'sec';
              'Mouth of Narr. Bay', 4,          'sec';
              'Upper East Passage', 17,         'seg';
              'Mid. East Passage',  18,         'seg';
              'Lower East Passage', 20, 'seg';
              'Upper West Passage', 66, 'seg';
              'Mid. West Passage',  65, 'seg';
              'Lower West Passage', 72,	'seg';
             };
%=========================================================================%

% Read shape file
shp_secs = shaperead(baysecs_shapefile);
shp_segs = shaperead(baysegs_shapefile);
grd_osom = grid_get(osom_gridfile);

% Determine grid points in each section and segment
n_sec    = size(sec_to_get,1);
sec_info = struct;
for i=1:n_sec
    
    % Section name
    sec_info(i).Name = sec_to_get{i,1};

    % Gather section shapes
    eval(['sec_info(i).Shape = shp_' sec_to_get{i,3} 's([' num2str(sec_to_get{i,2}) ']);']);

    % Polygon
    sec_info(i).Poly = simplify(polyshape( sec_info(i).Shape.X , sec_info(i).Shape.Y ));
    
end
clear i;

% Determine OSOM points within each polygon
isec = NaN(size(grd_osom.lon_rho));
for i=1:n_sec
    isec(grd_osom.mask_rho==1) = isinterior(sec_info(i).Poly,[grd_osom.lon_rho(grd_osom.mask_rho==1),grd_osom.lat_rho(grd_osom.mask_rho==1)]);
end

% Test
surf(grd_osom.lon_rho,grd_osom.lat_rho,isec,'edgecolor','none'); view(2);

% Save
save('osom_isec_nb.mat','isec','sec_info');