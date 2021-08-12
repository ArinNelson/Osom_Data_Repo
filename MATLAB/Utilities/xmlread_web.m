function xml_struct = xmlread_web(url_xml)
% Read xml-formatted URL from the web and output as a structure

    % Webread
    url_raw = webread(url_xml,weboptions('ContentType','text','Timeout',500));
    
    % Save as xml
    tmp_file = ['tmp_' sprintf('%12.0f',round(rand*1e12)) '.xml'];
    fid = fopen(tmp_file,'wt');
    fprintf(fid,url_raw);
    fclose(fid);
    
    % Read as xml
    xml_struct = xml2struct(tmp_file);
    
    % Delete temporary file
    delete(tmp_file);

end




%end