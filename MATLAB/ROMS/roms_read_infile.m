%function in = roms_read_infile(infile)
%=========================================================================%
% in = roms_read_infile(infile)
% reads in the variable names and values from the roms.in file specified by
% infile and returns them in a structure
% by Arin Nelson
% on 08/11/2021
% 
% NOTE: For now, this assumes only 1 grid is being used.  This code will
% eventually be updated to accomodate multiple/nested grids.
% 
% Last updated by Arin Nelson on 08/11/2021
%=========================================================================%

% FOR DEBUGGING
infile = 'D:\ROMS\Versions\DaveUllman\src\User\External\ocean.in';

% Open roms.in file
fid = fopen(infile,'rt');
if(fid==-1);    error(['Cannot find roms.in file: ' infile]);   end

% Read in lines of roms.in file and close it
inline = {};
while(~feof(fid));  inline{end+1} = fgetl(fid);     end
fclose(fid); clear fid;

% Remove comment lines (start with !) and blank lines
i = 1;
while(i<=numel(inline))
  if(isempty(inline{i}))
    inline(i) = [];
  else
    if( strcmp(inline{i}(1),'!') | all(isspace(inline{i})) )
      inline(i) = [];
    else
      i = i + 1;
    end
  end
end
clear i;

% Go through variables and get their values.
in= struct;
for i=1:numel(inline)

  % Split via spaces
  split_line = strsplit(inline{i},' ');
  
  % Remove empty entries
  n = 1;
  while(n<=numel(split_line))
    if(isempty(split_line{n}))
      split_line(n) = [];
    else
      n = n + 1;
    end
  end
  clear n;
  
  % If second entry is = or ==, is a variable
  % But don't save variables with ==
  if(numel(split_line)>2)
  if( (strcmp(split_line{2},'=') || strcmp(split_line{2},'==')) && ~contains(split_line{1},{'(',')'}) )
      
      % Variable value(s)
      tmp = cell(numel(split_line)-2,1);
      for j=1:numel(split_line)-2
          
          % If all characters, save as a string
          if(strcmp(split_line{j+2},'!'))
              tmp = tmp(1:j-1); break;
          else
          if(all(ismember(split_line{j+2}, '0123456789+-.*ed')))
              eval(['tmp{j} = ' split_line{j+2} ';']);
          else
              if(strcmp(split_line{j+2},'T'))
                  tmp{j} = true;
              elseif(strcmp(split_line{j+2},'F'))
                  tmp{j} = false;
              else
                tmp{j} = split_line{j+2};
              end
          end
          end
          
      end
      clear j;
      
      if(numel(tmp)==1)
          tmp = tmp{1};
      end
      
      % Save
      eval(['in.' split_line{1} '= tmp;']);
      clear tmp;
      
  end
  end

end
clear i;




%end