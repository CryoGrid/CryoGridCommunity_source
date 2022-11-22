%========================================================================
% CryoGrid FORCING class FORCING_tair_precip_mat
%
% simple model forcing for GROUND classes depending only on air temperature
% and precipitation (rain or snow).
%
% The data is obtained using the READ_FORCING_mat class. See this class for
% instructions about mat-file data format.
%
% The mandatory forcing variables are:
%
% Tair:      Air temperature (in degree Celsius)
% rainfall:  Rainfall (in mm/day)
% snowfall:  Snowfall (in mm/day)
%
% All forcing variables must be discretized identically, and one array of
% timestamps must be provided (t_stamp, in Matlab time / increment 1 
% corresponds to one day). 
%
% IMPORTANT POINT: the time series must be equally spaced in time, and this must be 
% really exact. When reading the timestamps from an existing data set (e.g. an Excel file),
% rounding errors can result in small differences in the forcing timestep, often less 
% than a second off. In this case, it is better to manually compile a new, equally spaced 
% timestep in Matlab.
%
% Authors
% T. Ingeman-Nielsen, October 2022
%========================================================================

classdef FORCING_tair_precip_mat < FORCING_base & READ_FORCING_mat
    
    
    methods

        function forcing = provide_PARA(forcing)         
            % INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.  

            forcing.PARA.filename = [];       % filename of Matlab file containing forcing data
			forcing.PARA.forcing_path = [];   % location (path) of forcing files
            forcing.PARA.start_time = [];     % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];       % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  % rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  % snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.all_rain_T = [];     % Temperature above which all precipitation is considered as rain
            forcing.PARA.all_snow_T = [];     % Temperature below which all precipitation is considered as snow
            forcing.PARA.heatFlux_lb = [];    % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];    % height above ground at which air temperature (and wind speed!) from the forcing data are applied.
            
        end
        
        function forcing = provide_CONST(forcing)
            
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
        function forcing = initialize_excel(forcing)
            
        end
 
        
        function forcing = finalize_init(forcing, tile)
            % FINALIZE_INIT  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
            
            variables = {'Tair'; 'rainfall'; 'snowfall'};
            [data, times] = read_mat([forcing.PARA.forcing_path forcing.PARA.filename], variables);
            
            for i=1:size(variables,1)
                if isfield(data, variables{i,1})
                    forcing.DATA.(variables{i,1}) = data.(variables{i,1});
                end
            end

            forcing.DATA.rainfall = data.rainfall.*forcing.PARA.rain_fraction;
            forcing.DATA.snowfall = data.snowfall.*forcing.PARA.snow_fraction;
            forcing.DATA.timeForcing = times;

            forcing = check_and_correct(forcing); % Remove known errors
            forcing = set_start_and_end_time(forcing); % assign start/end time
            
            %initialize TEMP
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.Tair=0;
            
        end


        function forcing = interpolate_forcing(forcing, tile)
            % Interpolate forcing data to timestep tile.t
            forcing = interpolate_forcing@FORCING_base(forcing, tile);

            forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
            forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
        end


        function fig = plot(forcing)
            TT = timetable(datetime(forcing.DATA.timeForcing,'ConvertFrom','datenum'), ...
                           forcing.DATA.Tair, ...
                           forcing.DATA.rainfall, ... 
                           forcing.DATA.snowfall, ...
                           'VariableNames', {'Tair', 'rainfall', 'snowfall'});
            TT.Properties.VariableUnits = {'degC', 'mm', 'mm'};
            stackedplot(TT);
            %hAx=gca;
            %hAx.TickLabelFormat='mm-yyyy';
            %datetick('x','mm-yyyy');
        end
        
        %-------------param file generation-----
        function forcing = param_file_info(forcing)
            forcing = provide_PARA(forcing);

            forcing.PARA.STATVAR = [];
            forcing.PARA.class_category = 'FORCING';
            
            forcing.PARA.comment.filename = {'filename of Matlab file containing forcing data'};
            
            forcing.PARA.default_value.forcing_path = {'forcing/'};
            forcing.PARA.comment.forcing_path = {'path where forcing data file is located'};
            
            forcing.PARA.comment.start_time = {'start time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.start_time.name =  'H_LIST';
            forcing.PARA.options.start_time.entries_x = {'year' 'month' 'day'};
            
            forcing.PARA.comment.end_time = {'end_time time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.end_time.name =  'H_LIST'; % 
            forcing.PARA.options.end_time.entries_x = {'year' 'month' 'day'};
            
            forcing.PARA.default_value.rain_fraction = {1};  
            forcing.PARA.comment.rain_fraction = {'rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)'};
            
            forcing.PARA.default_value.snow_fraction = {1};  
            forcing.PARA.comment.snow_fraction = {'snowfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)'};

            forcing.PARA.default_value.heatFlux_lb = {0.05};
            forcing.PARA.comment.heatFlux_lb = {'heat flux at the lower boundary [W/m2] - positive values correspond to energy gain'};
            
            forcing.PARA.default_value.airT_height = {2};  
            forcing.PARA.comment.airT_height = {'height above ground surface where air temperature from forcing data is applied'};

        end

    end

end