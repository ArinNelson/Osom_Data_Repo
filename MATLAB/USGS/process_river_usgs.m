function varargout = process_river_usgs(usgs_id,data_dir)
%=========================================================================%
% process_river_usgs(usgs_id, data_dir)
% Process river gauge data collected using 'gather_river_usgs.m' from 
% specified USGS river gauge ID usgs_id at the base data directory data_dir
% 
% stats = process_river_usgs(...0
% Returns the annual mean and variance of each available variable
% 
% by Arin Nelson
% on 08/23/2021
% 
% Last updated by Arin Nelson on 08/23/2021
%=========================================================================%

    % Check input args
    if(nargout>1); error('Too many output arguments.'); end

    % Variable info
    cfs_to_cms = 0.028316846592;    % convert cubic feet per second to cubic meters per second
    ft_to_m    = 1/3.281;           % convert feet to meters
    var_info   = {'00060',  'q',   cfs_to_cms; ...     % m3/s
                  '00065',  'h',   ft_to_m;    ...     % m
                  '00010',  'T',   1;          ...     % C
                  '00095',  'C',   1;          ...     % ?
                  '00300',  'DO2', 1;          ...     % mg/l
                  '00400',  'pH',  1;          ...     % ?
                 };
             
    % Load station info
	load([data_dir '/' usgs_id '/' usgs_id '_info.mat'],...
         'station_name','longitude','latitude','altitude','drainage_area');
     
    % Read list of available years
    tmp = dir([data_dir '/' usgs_id]);
    tmp = tmp([tmp.isdir]==1);
    tmp = tmp(3:end);
    tmp = cellfun(@str2num,{tmp.name});
    yy  = [tmp(:)];
	
    % Loop through years
    %yy = year_range(1) : 1 : year_range(end);
    ny = numel(yy);
    for iy=1:ny
        
        % List of available variables  
        file_list = ls([data_dir usgs_id '\' num2str(yy(iy)) '\*.mat']);
        if(~isempty(file_list)) 
        
            % Loop through variables
            for iv=1:size(file_list,1)
          
                % Variable name
                tmp = strsplit(file_list(iv,:),{'_','.'});  tmp=tmp{2};
                ii = find( strcmp(var_info(:,1),tmp)==1, 1 );
                if(~isempty(ii))
            
                    % Load data
                    tmp = load([data_dir usgs_id '\' num2str(yy(iy)) '\' file_list(iv,:)]);
          
                    % If variable doesn't yet exist, create it
                    var_name = var_info{ii,2};
                    if( ~exist(var_name,'var') )
                        eval([var_name '=[];']);
                        eval([var_name '_t=[];']);
                    end
                    
                    % Save data
                    nn = numel(tmp.time);
                    eval([var_name '(  end+1:end+nn) = tmp.value.*var_info{ii,3};']);
                    eval([var_name '_t(end+1:end+nn) = tmp.time;']);
                    clear nn;
                    
                end
                clear ii tmp;
                
            end
            clear iv;
            
        end
        clear file_list;
        
    end
    clear iy;
    
    % Quality control: gauge height
    % Outliers in h always above 10m (visually verified)
	if(exist('h','var'))
        h_t(h>10) = [];
        h(  h>10) = [];
    end
    
    % Quality control: flow rate
    % Low flow sometimes occur when gauge gets full of gunk
    % Also, missing data are 0's
	if(exist('q','var')) 
        %rvrgg(ig).q_t(rvrgg(ig).q<1) = [];
        %rvrgg(ig).q(  rvrgg(ig).q<1) = [];
        q_t(q==0) = [];
        q(  q==0) = [];
    end
    
    % Quality control: temperature
    if(exist('T','var'))
        
        % Sometimes units are *100
        if(max(T)>100)
            T = T ./ 100;
        end
            
        % If river temp <1, means sensor is covered in ice?
        T_t(T<1) = [];
        T(  T<1) = []; 
        
    end
    
    % Quality control, conductivity
    if(exist('C','var'))
       
        % none atm...
        
    end
    
    % Quality control, dissolved O2
    % High O2 usually means sensor is too close to surface
    if(exist('DO2','var'))
        DO2_t(DO2>60) = [];
        DO2(  DO2>60) = [];
    end
    
    % Save to file
    save_file = [data_dir '/' usgs_id '/' usgs_id '_master.mat'];
    save_str = 'save(save_file';
    if(exist('h','var'));   save_str = [save_str, ',''h'',''h_t'''];    end
    if(exist('q','var'));   save_str = [save_str, ',''q'',''q_t'''];    end
    if(exist('C','var'));   save_str = [save_str, ',''C'',''C_t'''];    end
    if(exist('T','var'));   save_str = [save_str, ',''T'',''T_t'''];    end
    if(exist('DO2','var')); save_str = [save_str, ',''DO2'',''DO2_t'''];    end
    save_str = [save_str ');'];
    eval(save_str);
    
    % compute stats if wanted
    if(nargout==1)
        var_name = {'q','h','T','C','DO2'};
        stats = zeros(ny,5,4);
        for iy=1:ny
        for iv=1:numel(var_name)
        if(exist(var_name{iv},'var'))
           eval(['t = ' var_name{iv} '_t;']);
           ii = find(year(t)==yy(iy)); 
           if(~isempty(ii))
               stats(iy,iv,1) = yy(iy);
               stats(iy,iv,2) = numel(ii);
               stats(iy,iv,3) = eval(['nanmean(' var_name{iv} '(ii));']);
               stats(iy,iv,4) = eval(['nanvar( ' var_name{iv} '(ii));']);
           end
        end
        end
        end
        varargout{1} = stats;
    end

%=========================================================================%