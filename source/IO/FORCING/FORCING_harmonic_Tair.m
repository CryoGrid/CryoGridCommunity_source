%========================================================================
% CryoGrid FORCING class FORCING_harmonic_tair
% simple model forcing for GROUND classes depending only on air temperature
% and precipitation (rain or snow).
% The class generates a harmonically oscillating air temperature
% timeseries, with a constant precipitation.
%
% The mandatory userdefined parameters are:
% start_time
% end_time
% time_step
% maat       Mean Annual Airt Temperature in degrees Celsius.
% amplitude  Amplitude of the yearly variation (degrees Celsius).
% lag        Phase lag (number of days to delay the harmonic oscillation). 
% period     The number of days in a year, defaults to 365.242 days.
% gradient   Gradual change in MAAT given in degrees/year, defaults to 0.
% precip:    Rain/snowfall in mm/day. Will be constant throughout the
%                simulation. Falls as snow when Tair<0, otherwise rain.
%
% T. Ingeman-Nielsen, December 2021
%========================================================================


classdef FORCING_harmonic_Tair < matlab.mixin.Copyable
    
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

            forcing.PARA.start_time = [];     % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];       % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.time_step = [];      % time step of timeseries [days] (e.g. 0.5 days)
            forcing.PARA.maat = [];           % Mean Annual Airt Temperature
            forcing.PARA.amplitude = [];      % Amplitude of the yearly variation
            forcing.PARA.lag = [];            % Phase lag (number of days to delay the harmonic oscillation)
            forcing.PARA.period = [];         % The number of days in a year, defaults to 365.242 days
            forcing.PARA.gradient = [];       % Gradual change in MAAT given in degrees/year
            forcing.PARA.precip = [];         % daily mean precipitation in mm/day
            forcing.PARA.heatFlux_lb = [];    % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];    % height above ground at which air temperature (and wind speed!) from the forcing data are applied.
        end
        
        function forcing = provide_CONST(forcing)
            
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
        function forcing = finalize_init(forcing, tile)
            % FINALIZE_INIT  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
            
            forcing.STATUS = 1;

            % handle start time
            forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            
            % handle end time
            forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));

            % generate time sequence
            forcing.DATA.timeForcing = forcing.PARA.start_time:forcing.PARA.time_step:forcing.PARA.end_time;
            forcing.PARA.time_zero = datenum(datetime(year(datetime(forcing.PARA.start_time,'ConvertFrom','datenum')),1,1,0,0,0));

            % handle period, if not specified
            if isempty(forcing.PARA.period) || isnan(forcing.PARA.period)
                forcing.PARA.period = 365.242;  % days in year
            end

            %initialize TEMP
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.Tair=0;
            
        end


        function forcing = interpolate_forcing(forcing, tile)
            % Interpolate forcing data to timestep tile.t
            t = tile.t;

            forcing.TEMP.Tair = forcing.PARA.maat - forcing.PARA.amplitude * cos(2 * pi * (t - forcing.PARA.time_zero - forcing.PARA.lag) / (forcing.PARA.period ));
            forcing.TEMP.Tair = forcing.TEMP.Tair + forcing.PARA.gradient * (t - forcing.PARA.start_time) / forcing.PARA.period;
            
            if forcing.TEMP.Tair <= 0
                forcing.TEMP.rainfall = 0;
                forcing.TEMP.snowfall = forcing.PARA.precip;
            else
                forcing.TEMP.snowfall = 0;
                forcing.TEMP.rainfall = forcing.PARA.precip;
            end
            
            forcing.TEMP.t = t;
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

    end    

end