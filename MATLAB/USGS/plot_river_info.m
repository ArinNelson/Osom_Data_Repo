%=========================================================================%
% plot_river_information.m
% Show river gauge locations and drainage basin(s) that they cover
% by Arin Nelson
% on 08/02/2021
% 
% Last modified by Arin Nelson on 08/02/2021
%=========================================================================%
clc; close all;

% Switches
Switch = zeros(9,1);
Switch(1) = 1;      % Load data
Switch(2) = 1;      % Plot

% Options
xLim = [-72 -70.7];
yLim = [41.3 42.36];
oceanColor = [0.3010 0.7450 0.9330];
landColor  = 'w';

% River information 
%              Name                  ID #       LON             LAT             DA      BA      RVR #
rvrgg_info = {'Taunton',            '01108000',	-70.9564307,	41.9339903,     261,	565;    ...
              'Blackstone',         '01113895',	-71.38144467,	41.8884331,     474,	540;    ...
              'Pawtuxet',           '01116500',	-71.4450575,	41.75093399,	200,	231.6;	...
              'Ten Mile',           '01109403',	-71.3503315,	41.83093377,	53.1,	54;     ...
              'Woonasquatucket',	'01114500',	-71.4872823,	41.8589884,     38.3,	50;     ...
              'Moshassuck',         '01114000',	-71.41061208,	41.833989,      23.1,	23.6;   ...
              'Palmer',             '01109220',	-71.2778279,	41.80926767,	30.9,	48.5;   ...
              'Hunt',               '01117000',	-71.4445016,	41.6412122,     22.9,	22.9;   ...
              'Mill',               '01108410', -71.09,         41.89972,       52.6,   NaN;    ...
              'Threemile',          '01109060',	-71.12282369,	41.8662122,     84.3,	NaN;    ...
              'Pawcatuck',          '01118500',	-71.8331247,	41.38371065,	295,    NaN;    ...
             };

% What to add to shoreline shape
shore_rvr = {'Ten Mile River'};
shore_lak = {'Blackstone River','Main Stem Pawtuxet River','Turner Reservoir (North & South)','Central Pond','Slater Park Pond'};
         
% Watershed info
ib   = [1 2 4];                         nb = numel(ib);     % Basins shapefile
it   = [5 6 7 10];                      nt = numel(it);     % HUC10 shapefile
iw   = [16 18 26 33 52 61 64 65 67];	nw = numel(iw);     % HUC12 shapefile
         
% GIS information
gis_dir        = 'F:/OSOM_Data_Repo/GIS/';
gis_states     = 'GeoData/States/statep010.shp';
gis_sections   = 'RIGIS/NBEP_Bay_Section_Boundaries/BAYSECTIONS_NBEP2017.shp';
gis_segments   = 'RIGIS/NBEP_Bay_Segment_Boundaries/BAYSEGMENTS_NBEP2017.shp';
gis_basins     = 'RIGIS/NBEP_River_Basin_Boundaries/BASINS_NBEP2017.shp';
gis_watersheds = 'RIGIS/NBEP_Watershed_Boundaries/Watershed_Boundary_Dataset%3A_HUC_12.shp';
gis_huc10      = 'RIGIS/NBEP_HUC10_Boundaries/HUC10_NBEP2017.shp';
gis_huc12      = 'RIGIS/NBEP_HUC12_Boundaries/HUC12_NBEP2017.shp';
gis_rivers     = 'RIGIS/RIGIS_Rivers_and_Streams/Rivers_and_Streams__RI_Integrated_Water_Quality_Monitoring_and_Assessment_Report_2012.shp';
gis_lakes      = 'RIGIS/RIGIS_Lakes_and_Ponds/Lakes_and_Ponds.shp';
gis_shorelines = 'RIGIS/RIGIS_CUSP/Continually_Updated_Shoreline_Product.shp';

%=========================================================================%
if(Switch(1))
    
    % Constants
    n_rvrgg = size(rvrgg_info,1);
    
    % Read in shape files
%     shp_states     = shaperead([gis_dir gis_states    ]);
%     shp_sections   = shaperead([gis_dir gis_sections  ]);
%     shp_segments   = shaperead([gis_dir gis_segments  ]);
    shp_basins     = shaperead([gis_dir gis_basins    ]);
%     shp_watersheds = shaperead([gis_dir gis_watersheds]);
    shp_huc10      = shaperead([gis_dir gis_huc10     ]);
    shp_huc12      = shaperead([gis_dir gis_huc12     ]);
%     shp_rivers     = shaperead([gis_dir gis_rivers    ]);
%     shp_lakes      = shaperead([gis_dir gis_lakes     ]);
%     shp_shorelines = shaperead([gis_dir gis_shorelines]);
   
%     % Rivers in shp_rivers
%     rvr_ndx = cell(n_rvrgg,1);
%     lak_ndx = cell(n_rvrgg,1);
%     for i=1:n_rvrgg
%         rvr_ndx{i} = find( contains({shp_rivers.NAME},rvrgg_info{i,1})==1 );
%         lak_ndx{i} = find( contains({shp_lakes.NAME },rvrgg_info{i,1})==1 );
%         if(~isempty(rvrgg_info{i,7}))
%         for j=1:numel(rvrgg_info{i,7})
%         	lak_ndx{i} = [lak_ndx{i}(:); find( contains({shp_lakes.NAME },rvrgg_info{i,7}{j})==1 )];
%         end
%         end
%     end
    
%     % Master RI shoreline map
%     shoreX = [shp_shorelines.X];
%     shoreY = [shp_shorelines.Y];
%     for i=1:numel(shore_rvr)
%         ii = find( strcmp( {shp_rivers.NAME}, shore_rvr{i})==1 );
%         shoreX = [shoreX(:)', shp_rivers(ii).X];
%         shoreY = [shoreY(:)', shp_rivers(ii).Y];
%     end
%     for i=1:numel(shore_lak)
%         ii = find( strcmp( {shp_lakes.NAME}, shore_lak{i})==1 );
%         shoreX = [shoreX(:)', shp_lakes(ii).X];
%         shoreY = [shoreY(:)', shp_lakes(ii).Y];
%     end
%     

% NOTE: THESE INDICES DON'T SEEM TO MATCH THE ONES ONLINE...
% 
%     lak_id = [82, 2517, 2163, 2164, ... % Blackstone
%               1228, 1262, 1269, 1599, 1651, 1800, 3735, ... % Ten Mile
%               360, 462, 479, 513, 524, 526, 1115, 1243, 1251, 1260, 1997, ... % Moshassuck
%               377, 463, 514, 617, 695, 710, 1642, 2004, 2005, ... % Woonasquatucket
%               1844, 1869, 1877, ... % Main Pawtuxet
%               665, 942, 1989, 1725, 1808, ... % North Pawtuxet
%               1151, 1158, 1165, 1411, 2159, 2203, 2572, 2583, 2611, 2693, 1885, 2746, 2782, 2158, ... % South Pawtuxet
%               2909, ... % Hunt
%               2813, 1282, 2796, ... % Maskerchugg
%               1161, 1169, 2678, 2683, 2687, ... % Hardig
%               4, 2160, 2161, 2162, 2678, 2683, 2687, 3401, 3410, 3427, 3519, 3533]; % Pawcatuck
%      rvr_id = [3371, 3377, 3378, 3382, 3844, 3399, 3412, 3418, 3437, 3440, 3443, 3461, 3462, 3463, ... % Hunt
%                 ];         
    

% Basin and HUC10/12 ID's
basin_id = [1, 2, 4, 5]; 
huc10_id = [5, 6, 7,  10, 12]; 
huc12_id = [5, 7, 16, 18, 26, 33, 39, 50, 64, 65]; 

% Plotting info
cmap = lines(11);
basin_info = {'Blackstone R.',      'basins', 1,'cmap(1,:)'; ...
              'Pawcatuck R.',       'basins', 2,'cmap(2,:)'; ...
              'Pawtuxet R.',        'basins', 4,'cmap(4,:)'; ...
              'Coastal Ponds',      'basins', 5,'[0.75 0.75 0.75]'; ...
              'Taunton R. (N)',     'huc10',  5,'cmap(3,:)'; ...
              'Taunton R. (S)',     'huc10',  6,'cmap(3,:)'; ...
              'Threemile R.',       'huc10',  7,'cmap(6,:)'; ...
              '(riparian)',         'huc10', 12, '[0.75 0.75 0.75]'; ...
              'Ten Mile R.',        'huc10', 10,'cmap(5,:)'; ...
              'Assonet R.',         'huc12',  5,'[0.75 0.75 0.75]'; ...
              'Greenwich Bay',      'huc12', 16,'cmap(8,:)'; ...
              'Hunt R.',            'huc12', 18,'cmap(8,:)'; ...
              'Moshassuck R.',      'huc12', 26,'cmap(7,:)'; ...
              'Palmer R.',          'huc12', 33,'cmap(9,:)'; ...
              'Warren R.',          'huc12',  7,'cmap(9,:)'; ...
              'Quequechan R.',      'huc12', 39,'[0.75 0.75 0.75]'; ...
              '(riparian)',         'huc12', 50,'[0.75 0.75 0.75]'; ...
              'Woonasquatucket R.', 'huc12', 64,'cmap(10,:)'; ...
              'Mill R.',            'huc12', 65,'cmap(11,:)'; ...
             };

end
%=========================================================================%
if(Switch(2))
    
    % Figure
    fig = figure('units','normalized','outerposition',[0 0 1 1],'defaultaxesfontsize',12);
    
    % Axis
    %axx = usamap({'ct','ma','ri','ny'});
    %setm(axx, 'FFaceColor', oceanColor)
    axx = axes('units','normalized','position',[0 0.075 1 0.90]);
    hold(axx,'on');
    
    % Show states
    %plot(shoreX,shoreY,'-b');
    %ii_ri = find( strcmp( {shp_states.STATE}, 'Rhode Island'  )==1 );   mapshow(shp_states(ii_ri),'FaceColor',landColor);
    %plot([shp_shorelines.X],[shp_shorelines.Y],'-k');
    %ii_ma = find( strcmp( {shp_states.STATE}, 'Massachusetts' )==1 );   mapshow(shp_states(ii_ma),'FaceColor',landColor);
    %ii_ct = find( strcmp( {shp_states.STATE}, 'Connecticut'   )==1 );	mapshow(shp_states(ii_ct),'FaceColor',landColor);
    %ii_ny = find( strcmp( {shp_states.STATE}, 'New York'      )==1 );   mapshow(shp_states(ii_ny),'FaceColor',landColor); 
   
    % Plot major rivers
    %plot([shp_rivers(rvr_id).X],[shp_rivers(rvr_id).Y],'-b'); 
    %plot([shp_lakes( lak_id).X],[shp_lakes( lak_id).Y],'-b'); 
%     mapshow(shp_rivers(rvr_id));
%     mapshow(shp_lakes( lak_id));
    
%     % Show major river basins
%     plot([shp_huc10(huc10_id).X],[shp_huc10(huc10_id).Y],'-k');
%     plot([shp_huc12(huc12_id).X],[shp_huc12(huc12_id).Y],'-k');
%     plot([shp_basins(basin_id).X],[shp_basins(basin_id).Y],'-k');
    for i=1:size(basin_info,1)
       eval(['mapshow(shp_' basin_info{i,2} '(' num2str(basin_info{i,3}) '),''FaceColor'',' basin_info{i,4} ');']);
       %eval(['xx=prctile(shp_' basin_info{i,2} '(' num2str(basin_info{i,3}) ').X,50);']);
       %eval(['yy=prctile(shp_' basin_info{i,2} '(' num2str(basin_info{i,3}) ').Y,50);']);
       %text(xx,yy,basin_info{i,1});
    end
   
    % Show river gauge locations
    for i=1:n_rvrgg
        plot(rvrgg_info{i,3},rvrgg_info{i,4},'kp','markersize',16,'color','k','MarkerFaceColor','c');
    end
    
    % Plot properties
    hold(axx,'off'); 
    axis equal; axis tight; box on;
    xlim(xLim); ylim(yLim); xlabel('Longitude (^oE)'); ylabel('Latitude (^oN)');
    title('Rhode Island & Narragansett Bay River Basins & Gauges');
    %set(gca,'color',oceanColor);
    
    % Legend using invisible axis
    
end