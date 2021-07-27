function [var_name, var_value] = webread_usgs(str_url)

  % Read data from url (assumedly in rdb format)
  txt_raw = webread(str_url,weboptions('Timeout', 120));
  
  % Split at newline characters
  txt_lines = strsplit(txt_raw,'\n');
  
  % Remove empty lines and header lines
  i = 1;
  while i<=numel(txt_lines)
    if(isempty(txt_lines{i}))
      txt_lines(i) = [];
    else
      if(txt_lines{i}(1)=='#')
        txt_lines(i) = [];
      else
        i = i + 1;
      end
    end
  end
  clear i;
  
  % First line is variable names
  var_name     = strsplit( txt_lines{1}, '\t' );
  txt_lines(1) = [];
  
  % Second line is time to retrieve, so it can be ignored
  txt_lines(1) = [];
  
  % Remaining lines are data values
  if(numel(txt_lines)>0)
    var_value = cell(numel(txt_lines),numel(var_name));
    for i=1:numel(txt_lines)
    
      % Split this line  
      this_line = strsplit(txt_lines{i},'\t');
    
      % Save the vaues
      if(~all(cellfun(@isempty,this_line)))
      for j=1:numel(this_line)
      if(~isnan(str2double(this_line{j})))    
        var_value{i,j} = str2double(this_line{j});
      else
        var_value{i,j} = this_line{j};
      end
      end
      end
    
    end
    clear i j this_line;
  else
    var_value = '';
  end
  
  % Remove constant values
  ii = find( strcmp( var_name, 'agency_cd' )==1 );    var_name(ii) = [];  var_value(:,ii) = [];     clear ii;
  ii = find( strcmp( var_name, 'site_no'   )==1 );    var_name(ii) = [];  var_value(:,ii) = [];     clear ii;
  ii = find( contains( var_name, '_cd' )==1     );    var_name(ii) = [];  var_value(:,ii) = [];     clear ii;
  
end