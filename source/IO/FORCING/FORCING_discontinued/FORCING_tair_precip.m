%========================================================================
% CryoGrid FORCING class FORCING_tair_precip
% simple model forcing for GROUND classes depending only on air temperature
% and precipitation (rain or snow).
% The data must be stored in a Matlab “.mat” file which contains 
% a struct FORCING with field “data”, which contain the time series of the actual 
% forcing data, e.g. FORCING.data.Tair contains the time series of air temperatures. 
% Have a look at the existing forcing files in the folder “forcing” and prepare 
% new forcing files in the same way. 
%
% The mandatory forcing variables are:
% Tair:      Air temperature (Tair, in degree Celsius)
% rainfall:  Rainfall (rainfall, in mm/day), 
% snowfall:  Snowfall (snowfall, in mm/day)
% t_span:    Timestamp (t_span, in Matlab time / increment 1 corresponds to one day)
%
% IMPORTANT POINT: the time series must be equally spaced in time, and this must be 
% really exact. When reading the timestamps from an existing data set (e.g. an Excel file),
% rounding errors can result in small differences in the forcing timestep, often less 
% than a second off. In this case, it is better to manually compile a new, equally spaced 
% timestep in Matlab.
%
% T. Ingeman-Nielsen, S. Westermann, J. Scheer, December 2021
%========================================================================

classdef FORCING_tair_precip < matlab.mixin.Copyable
    
    properties
        forcing_index
        DATA            % forcing data time series
        TEMP            % forcing data interpolated to a timestep
        PARA            % parameters
        STATUS         
        CONST         
    end
    
    
    methods

        function forcing = provide_PARA(forcing)         
            % INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.  

            forcing.PARA.filename = [];       % filename of Matlab file containing forcing data
			forcing.PARA.forcing_path = [];   % location (path) of forcing files
            forcing.PARA.latitude = [];       % 
            forcing.PARA.longitude = [];      % 
            forcing.PARA.altitude = [];       % 
            forcing.PARA.start_time = [];     % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];       % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  % rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  % snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
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
            
			temp = load([forcing.PARA.forcing_path forcing.PARA.filename], 'FORCING');
            
            forcing.DATA.Tair = temp.FORCING.data.Tair;
            forcing.DATA.rainfall=temp.FORCING.data.rainfall.*forcing.PARA.rain_fraction;
            forcing.DATA.snowfall=temp.FORCING.data.snowfall.*forcing.PARA.snow_fraction;
            forcing.DATA.timeForcing = temp.FORCING.data.t_span;
            
            if std(forcing.DATA.timeForcing(2:end,1)-forcing.DATA.timeForcing(1:end-1,1))>=1e-10
                error('timestamp of forcing data is not in regular intervals -> check, fix and restart')
            else
                forcing.STATUS = 1;
            end

            % handle start time, if specified
            if isempty(forcing.PARA.start_time) || ~ischar(forcing.PARA.start_time)
                forcing.PARA.start_time = forcing.DATA.timeForcing(1,1);
            else
                forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            end
            
            % handle end time, if specified
            if isempty(forcing.PARA.end_time) || isnan(forcing.PARA.end_time(1,1))
                forcing.PARA.end_time = floor(forcing.DATA.timeForcing(end,1));
            else
                forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));
            end
            
            %initialize TEMP
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.Tair=0;
            
        end

        function forcing = interpolate_forcing(forcing, tile)
            % Interpolate forcing data to timestep tile.t
            t = tile.t;
            times = forcing.DATA.timeForcing;
            posit = floor((t-times(1,1))./(times(2,1)-times(1,1)))+1;

            forcing.TEMP.snowfall = forcing.lin_interp(t, posit, times, forcing.DATA.snowfall);
            forcing.TEMP.rainfall = forcing.lin_interp(t, posit, times, forcing.DATA.rainfall);
            forcing.TEMP.Tair =     forcing.lin_interp(t, posit, times, forcing.DATA.Tair);

            forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
            forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
            forcing.TEMP.t = t;
        end
% 
%         function xls_out = write_excel(forcing)
% 			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
% 			error('This function is not implemented/updated for this specific class')
%             xls_out = {'FORCING','index',NaN,NaN;'FORCING_seb',1,NaN,NaN;NaN,NaN,NaN,NaN;'filename',NaN,NaN,NaN;'start_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the first timestamp of the forcing data set will be used';'end_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the last timestamp of the forcing data set will be used';'rain_fraction',1,'[-]','rainfall in forcing file multiplied by this number';'snow_fraction',1,'[-]','snowfall in forcing file multiplied by this number';'latitude',NaN,'[degree]','geographical coordinates';'longitude',NaN,'[degree]',NaN;'altitude',NaN,'[m]','a.s.l.';'domain_depth',100,'[m]','should match a GRID point, model domain extends to this depth';'heatFlux_lb',0.0500000000000000,'[W/m2]','geothermal heat flux';'airT_height',2,'[m]','height of air temperature';'FORCING_END',NaN,NaN,NaN};
%         end

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

    
    methods(Static)
        function value = lin_interp(t, posit, times, data)
			% t       is the current time
			% times   is the vector of times at which forcing data are available
			% data    is the vector of forcing data (on parameter)
			% posit   is an index into the time vector to the largest specified time step before current time 
            value = data(posit,1) + (data(posit+1,1) - data(posit,1)).*(t-times(posit,1))./(times(2,1)-times(1,1));
        end        
    end
    

end