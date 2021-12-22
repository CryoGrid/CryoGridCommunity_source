%========================================================================
% CryoGrid FORCING class FORCING_ubT_lbT_fixed
% simple model forcing for GROUND classes depending only on fixed upper
% and lower boundary temperatures
%
% The mandatory userdefined parameters are:
% start_time
% end_time
% time_step
% T_ub   Upper boundary temperature.
% T_lb   Lower boundary temperature.
%
% T. Ingeman-Nielsen, December 2021
%========================================================================


classdef FORCING_ubT_lbT_fixed < matlab.mixin.Copyable
    
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
            forcing.PARA.T_ub = [];           % 
            forcing.PARA.T_lb = [];           % 
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
            
            forcing.STATUS = 1;

            % handle start time
            forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            
            % handle end time
            forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));

            % generate time sequence
            forcing.DATA.timeForcing = forcing.PARA.start_time:forcing.PARA.time_step:forcing.PARA.end_time;
            forcing.PARA.time_zero = datenum(datetime(year(datetime(forcing.PARA.start_time,'ConvertFrom','datenum')),1,1,0,0,0));

            %initialize TEMP
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.Tair=0;
            forcing.TEMP.T_ub=0;
            forcing.TEMP.T_lb=0;
            
        end


        function forcing = interpolate_forcing(forcing, tile)
            % Interpolate forcing data to timestep tile.t
            t = tile.t;
            
            forcing.TEMP.T_ub = forcing.PARA.T_ub;
            forcing.TEMP.T_lb = forcing.PARA.T_lb;
            forcing.TEMP.Tair = forcing.TEMP.T_ub;
            forcing.TEMP.t = t;
        end

        function xls_out = write_excel(forcing)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			error('This function is not implemented/updated for this specific class')
            xls_out = {'FORCING','index',NaN,NaN;'FORCING_seb',1,NaN,NaN;NaN,NaN,NaN,NaN;'filename',NaN,NaN,NaN;'start_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the first timestamp of the forcing data set will be used';'end_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the last timestamp of the forcing data set will be used';'rain_fraction',1,'[-]','rainfall in forcing file multiplied by this number';'snow_fraction',1,'[-]','snowfall in forcing file multiplied by this number';'latitude',NaN,'[degree]','geographical coordinates';'longitude',NaN,'[degree]',NaN;'altitude',NaN,'[m]','a.s.l.';'domain_depth',100,'[m]','should match a GRID point, model domain extends to this depth';'heatFlux_lb',0.0500000000000000,'[W/m2]','geothermal heat flux';'airT_height',2,'[m]','height of air temperature';'FORCING_END',NaN,NaN,NaN};
        end

        function fig = plot(forcing)
			error('This function is not implemented/updated for this specific class')
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