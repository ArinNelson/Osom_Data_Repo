function ggl_table = webread_googlesheet(ggl_id,sheet_name)
%--------------------------------------------------------------------------
% webread_googlesheet - read in a Google Sheet and save to a table
% by Arin Nelson
% last updated 09/30/2021
% 
% usage: ggl_table = webread_googlesheet(ggl_id, sheet_name)
% 
%   ggl_id is the ID number from the Google sheet's url
%   note: url format is https://docs.google.com/spreadsheets/d/<ggl_id>
% 
%   sheet_name is the name of the sheet to read in
% 
%--------------------------------------------------------------------------

    % URL of the Google Sheet
    ggl_url   = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',ggl_id,sheet_name);
    
    % Read it in as a table
    ggl_table = webread(ggl_url);
    
end