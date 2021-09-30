function data_table = webread_usgs(str_url)

    % To read USGS data into a table...
    readtableweb = @(filename)readtable(filename,'CommentStyle','#');
    myoptions    = weboptions('ContentReader',readtableweb,'Timeout',600);

    % Read raw tab-delimited text
    data_table = webread(str_url,myoptions);

    % First row is a timing string
    data_table = data_table(2:end,:);
    
end