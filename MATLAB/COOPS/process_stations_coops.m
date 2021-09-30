%=========================================================================%
% Process land/dock station data from CO-OPS
% Arin Nelson
% on 08/01/2021
% 
% Last edited 08/01/2021
%=========================================================================%
clc; clear mex; addpath('../Utilities');

% Options
coops_dir = 'F:/OSOM_Data_Repo/Stations/COOPS/';
coops_id  = {'8454049'};
for i=1:numel(coops_id)
   process_station_coops(coops_dir,coops_id{i}); 
end




