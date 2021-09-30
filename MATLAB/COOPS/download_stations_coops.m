%=========================================================================%
% Gather land/dock station data from CO-OPS
% Arin Nelson
% on 07/23/2021
% 
% Last edited 08/01/2021
%=========================================================================%
clear; clc; clear mex; addpath('../Utilities');

% Options
get_year = 1990:2020;
save_dir = 'F:/OSOM_Data_Repo/Stations/COOPS/';
sttn_id  = {'8452944','8454049','8452660','8447386'};

% Constants
var_to_get = {'water_level','water_temperature','salinity',...
              'air_temperature','air_pressure'}; % Variables to look for       

%=========================================================================%
for is=1:numel(sttn_id)
clc; disp(['On station ' num2str(is) '/' num2str(numel(sttn_id)) '...']);    
    download_station_coops([save_dir '/COOPS_' sttn_id{is} '.nc'],sttn_id{is},var_to_get);
end
%=========================================================================%