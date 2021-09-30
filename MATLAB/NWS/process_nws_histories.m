nws_dir   = 'F:\OSOM_Data_Repo\Stations\NWS\';
nws_input = ls([nws_dir '*_raw.nc']);
nws_output = nws_input(:,1:8);
for i=1:size(nws_output,1)
    process_nws_history([nws_dir nws_input(i,:)],[nws_dir nws_output(i,:) 'daily.nc']);
end
