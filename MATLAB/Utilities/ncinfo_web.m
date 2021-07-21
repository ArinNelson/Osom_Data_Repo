function [vardata]=ncinfo_web(max_time,ncFile)
% function [vardata]=ncinfo_web(max_time,ncFile)
% wrapper function for ncinfo to be used when reading netcdf files from URL
% (not local file). this function uses parallel processing capability to
% monitor the ncread request and to kill it (and retry) if taking too long.
%  1st input argument is maximum time (seconds) to wait for an ncread
%  request to complete. 
%  other input arguments are identical to those for
%  ncread.
% NOTE that this function will run indefinitely if the request is never
% fulfilled. probably should put in a maximum number of attempts.
% 
% if max_time is 0, no wrapper is used (i.e., same as calling ncinfo)
% 
% modified from ncread_web.m provided by David Ullman, URI
if(max_time==0)
  vardata = ncinfo(ncFile);
else
  got_data=0;
  var_data=[];
  while (~got_data)
	F=parfeval(@ncinfo,1,ncFile);
	tic
	while( toc <= max_time )
    if( strfind(F.State,'finished') )
	  vardata=fetchOutputs(F);
	  got_data=1;
	  break
    end
    end
	cancel(F)
    got_data=1;
  end
end

end
