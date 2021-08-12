%=========================================================================%
% process_rivers_usgs.m
% Gather data from multiple rivers
% 
% Arin Nelson
% on 07/29/2021
% 
% Last edited 07/29/2021
%=========================================================================%

% Options
rvrgg_year = [1990 2020];     % Collect data through these years
rvrgg_dir  = 'D:/OSOM_Data_Repo/USGS/RiverGauges/';
% sttn_dir   = 'D:/OSOM_Data_Repo/Stations/COOPS/';
% 
% % Google doc with meta data
% ggl_base   = ['https://docs.google.com/spreadsheets/d/'      ...
%               '1Uc7mBJH-hggEDFPP6eWMh6tITHHFbDG2jtcLnI0ubdo' ...
%               '/gviz/tq?tqx=out:csv&sheet='];
% ggl_rvrgg  = 'River Sources';
% ggl_sttn   = 'Stations';
          
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

%=========================================================================%

	% Stations available
    rvrgg_list = ls([rvrgg_dir '/0*']);
    n_rvrgg    = size(rvrgg_list,1);
    
    % Load info from each station
    rvrgg = struct;
    for i=1:n_rvrgg
        tmp = load([rvrgg_dir '/' rvrgg_list(i,:) '/' rvrgg_list(i,:) '_info.mat']);
        rvrgg(i).id   = rvrgg_list(i,:);
        rvrgg(i).name = tmp.station_name;
        rvrgg(i).lon  = tmp.longitude;
        rvrgg(i).lat  = tmp.latitude;
        rvrgg(i).alt  = tmp.altitude;
        rvrgg(i).area = tmp.drainage_area;
    end
    clear i tmp;
    
    % Lookit?
    %{rvrgg.id; rvrgg.name; rvrgg.lon; rvrgg.lat; rvrgg.area}'
    
    % Gather available data from each station during the provided year span
    yy = rvrgg_year(1) : rvrgg_year(end);
    ny = numel(yy);
    for ig=1:n_rvrgg
    for iy=1:ny
    if(exist( [rvrgg_dir rvrgg(ig).id '/' num2str(yy(iy))], 'dir' )==7)
        
      % List of available variables  
      file_list = ls([rvrgg_dir rvrgg(ig).id '\' num2str(yy(iy)) '\*.mat']);
      if(~isempty(file_list))  
      
        % Loop through variables
        for iv=1:size(file_list,1)
          
            % Variable name
            tmp = strsplit(file_list(iv,:),{'_','.'});  tmp=tmp{2};
            ii = find( strcmp(var_info(:,1),tmp)==1 );
            if(~isempty(ii))
      
                % Load data
                tmp = load([rvrgg_dir rvrgg(ig).id '\' num2str(yy(iy)) '\' file_list(iv,:)]);
          
                % If variable doesn't yet exist, create it
                var_name = var_info{ii,2};
                if( ~isfield(rvrgg(ig),var_name) )
                    eval(['rvrgg(ig).' var_name '=[];']);
                    eval(['rvrgg(ig).' var_name '_t=[];']);
                end
                
                % Save variables
                nn = numel(tmp.time);
                eval(['rvrgg(ig).' var_name '(  end+1:end+nn) = tmp.value.*var_info{ii,3};']);
                eval(['rvrgg(ig).' var_name '_t(end+1:end+nn) = tmp.time;']);
                
                
            end
            clear ii tmp;
            
        end
        clear iv;
        
      end
      clear file_list;
      
    end
    end
    end
    clear ig iy yy ny;
    
%-------------------------------------------------------------------------%
        
    % Perform quality control
    for ig=1:n_rvrgg
    
        % Outliers in h always above 10m (visually verified)
        if(~isempty(rvrgg(ig).h))
            rvrgg(ig).h_t(rvrgg(ig).h>10) = [];
            rvrgg(ig).h(  rvrgg(ig).h>10) = [];
        end
        
        % Missing values in q always 0
        if(~isempty(rvrgg(ig).q))
            rvrgg(ig).q_t(rvrgg(ig).q==0) = [];
            rvrgg(ig).q(  rvrgg(ig).q==0) = [];
        end
        
        % River stuff
        if(~isempty(rvrgg(ig).T))
            
            % Sometimes units are *100
            if(max(rvrgg(ig).T)>100)
                rvrgg(ig).T = rvrgg(ig).T ./ 100;
            end
            
            % If river temp <1, usually means its ice-covered
            rvrgg(ig).T_t(rvrgg(ig).T<1) = [];
            rvrgg(ig).T(  rvrgg(ig).T<1) = []; 
            
        end
        
        % Conductivity stuff
%         if(~isempty(rvrgg(ig).C))
%          
% %             % Sometimes units are *10
% %             if(max(rvrgg(ig).C)>500)
% %                 rvrgg(ig).C = rvrgg(ig).C./100;
% %             end
% %          
% %             % Convert conductivity to salinity
% %             % Looks like they have different units...
% %             %rvrgg(ig).C = cond_to_saln(rvrgg(ig).C,25);
%             
%         end

        % High O2 usually means sensor is too close to surface
        if(~isempty(rvrgg(ig).DO2))
            rvrgg(ig).DO2_t(rvrgg(ig).DO2>60) = [];
            rvrgg(ig).DO2(  rvrgg(ig).DO2>60) = [];
        end
        
    end
    clear ig;
    
    % Lookit?
    %hold on; for i=1:n_rvrgg; plot(rvrgg(i).DO2_t,rvrgg(i).DO2,'.'); end
            
%-------------------------------------------------------------------------%

    % Save river gauge data into a master file per station
    for ig=1:n_rvrgg
        
        % Save file
        save_file = [rvrgg_dir '/' rvrgg(ig).id '/' rvrgg(ig).id '_master.mat'];
        
        % Save variables
        rvrgg_flds = fieldnames(rvrgg(ig));
        save_str = 'save(save_file';
        for i=1:numel(rvrgg_flds)
          eval([rvrgg_flds{i} '=rvrgg(ig).' rvrgg_flds{i} ';']);
           save_str = [save_str ',''' rvrgg_flds{i} '''']; 
        end
        save_str = [save_str ');'];
       
        % Save
        eval(save_str);
        
        % Clean-up
        clear save_file rvrgg_flds save_str;

    end
    clear ig;

%-------------------------------------------------------------------------%