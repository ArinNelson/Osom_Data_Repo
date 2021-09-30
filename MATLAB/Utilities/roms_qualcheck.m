%=========================================================================%
% roms_qualcheck.m
% Quality check the specified input file by:
%   finding missing data
%   interpolating missing data in time
% 
% last edited by Arin Nelson on 07/28/2021
%=========================================================================%
clc; clear mex;

% fname = 'OSOM_bry_2018_DOPPIO_best.nc';
% info  = ncinfo(fname);
% for i=1:numel(info.Variables)
%     test = ncread(fname,info.Variables(i).Name);
%     if(any(test>1e36))
%         pause(1);
%     end
% end

% Truncate in time
it = 1:8670;
nt = numel(it);
info=ncinfo('OSOM_bry_2018_DOPPIO_best.nc');
info.Dimensions(end).Length = nt;
for i=1:numel(info.Variables)
    if(~isempty(info.Variables(i).Dimensions))
    if(any( strcmp({info.Variables(i).Dimensions.Name},'bry_time')==1 ))
        info.Variables(i).Dimensions(end).Length = nt;
        info.Variables(i).Size(end) = nt;
    end
    end
end
ncwriteschema('OSOM_bry_2018_DOPPIO_best_fix.nc',info);
for i=1:numel(info.Variables)
    if(~isempty(info.Variables(i).Dimensions))
    if(any( strcmp({info.Variables(i).Dimensions.Name},'bry_time')==1 ))
        tmp = ncread('OSOM_bry_2018_DOPPIO_best.nc',info.Variables(i).Name);
        str = '(';
        n   = 1;
        if(numel(size(tmp))==2 & any(size(tmp)==1))
            str = '(';
        else
          while n<numel(size(tmp))
            str = [str ':,'];
            n   = n + 1;
          end
        end
        str = [str 'it);'];
        eval(['tmp = tmp' str]);
        ncwrite('OSOM_bry_2018_DOPPIO_best_fix.nc',info.Variables(i).Name,tmp);
        clear tmp str n;
    else
        ncwrite('OSOM_bry_2018_DOPPIO_best_fix.nc',info.Variables(i).Name,ncread('OSOM_bry_2018_DOPPIO_best.nc',info.Variables(i).Name));
    end
    else
        ncwrite('OSOM_bry_2018_DOPPIO_best_fix.nc',info.Variables(i).Name,ncread('OSOM_bry_2018_DOPPIO_best.nc',info.Variables(i).Name));
    end
end








