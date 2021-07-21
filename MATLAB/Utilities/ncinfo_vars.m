function info_var = ncinfo_vars(ncfile)
% Parse out variable names, dimensions, and units from NetCDF

info=ncinfo(ncfile);
info_var = cell(numel(info.Variables),3);
for i=1:numel(info.Variables)
  info_var{i,1} = info.Variables(i).Name;
  info_var{i,2} = info.Variables(i).Size;
  ii = find( strcmp({info.Variables(i).Attributes.Name},'units')==1 );
  if(~isempty(ii))
    info_var{i,3} = info.Variables(i).Attributes(ii).Value;
  end
end
clear i ii;




end