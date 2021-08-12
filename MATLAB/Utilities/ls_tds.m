function file_list = ls_tds(url_catalog)
% Read list of files from THREDDS catalog at specified catalog.html

    % DEBUGGING
    url_catalog = 'https://www.ncei.noaa.gov/thredds/catalog/model-namanl-old/201001/20100101/catalog.xml';

    % Get xml structure from this url
    xml_catalog = xmlread_web(url_catalog);
    
    % Parse out file URLs
    url_list  = xml_catalog.catalog.dataset.dataset;
    n_url     = numel(url_list);
    file_list = cell(n_url,2);
    for i=1:n_url
       file_list{i,1} = url_list{i}.Attributes.urlPath;
       file_list{i,2} = url_list{i}.Attributes.name;
    end
    clear i;

end