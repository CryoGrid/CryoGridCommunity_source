%========================================================================
% CryoGrid FORCING class FORCING_Tair_Tlb_fixed
% simple model forcing for GROUND classes depending only on an upper
% boundary temperature (Tair) forcing, and a fixed lower boundary
% temperature.
% The forcing data must be stored in a Matlab “.mat” file which contains 
% a struct FORCING with field “DATA”, which contain the time series of the actual 
% forcing data, e.g. FORCING.data.Tair contains the time series of air temperatures. 
% Have a look at the existing forcing files in the folder “forcing” and prepare 
% new forcing files in the same way. 
%
% The mandatory forcing variables are:
% Tair:      Air temperature (Tair, in degree Celsius)
% t_span:    Timestamp (as Matlab datenum / increment 1 corresponds to one day)
%
% T_lb is specified in the parameter input file, not in the forcing mat file.
%
% T. Ingeman-Nielsen, S. Westermann, December 2021
%========================================================================

classdef FORCING_Tair_Tlb_fixed < matlab.mixin.Copyable
    
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
            forcing.PARA.start_time = [];     % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];       % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.T_lb = [];           % Fixed temperature at lower boundary
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
            forcing.TEMP.Tair = 0;
			forcing.TEMP.T_lb = forcing.PARA.T_lb;
            
        end

        function forcing = interpolate_forcing(forcing, tile)
            % Interpolate forcing data to timestep tile.t
            t = tile.t;
            times = forcing.DATA.timeForcing;
            posit = floor((t-times(1,1))./(times(2,1)-times(1,1)))+1;

            forcing.TEMP.Tair = forcing.lin_interp(t, posit, times, forcing.DATA.Tair);
			forcing.TEMP.T_lb = forcing.PARA.T_lb;
            forcing.TEMP.t = t;
        end

        function xls_out = write_excel(forcing)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			error('This function is not implemented/updated for this specific class')
            xls_out = {'FORCING','index',NaN,NaN;'FORCING_seb',1,NaN,NaN;NaN,NaN,NaN,NaN;'filename',NaN,NaN,NaN;'start_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the first timestamp of the forcing data set will be used';'end_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the last timestamp of the forcing data set will be used';'rain_fraction',1,'[-]','rainfall in forcing file multiplied by this number';'snow_fraction',1,'[-]','snowfall in forcing file multiplied by this number';'latitude',NaN,'[degree]','geographical coordinates';'longitude',NaN,'[degree]',NaN;'altitude',NaN,'[m]','a.s.l.';'domain_depth',100,'[m]','should match a GRID point, model domain extends to this depth';'heatFlux_lb',0.0500000000000000,'[W/m2]','geothermal heat flux';'airT_height',2,'[m]','height of air temperature';'FORCING_END',NaN,NaN,NaN};
        end

        function fig = plot(forcing)
            TT = timetable(datetime(forcing.DATA.timeForcing,'ConvertFrom','datenum'), ...
                           forcing.DATA.Tair, ...
                           'VariableNames', {'Tair'});
            TT.Properties.VariableUnits = {'degC'};
            stackedplot(TT);
            %hAx=gca;
            %hAx.TickLabelFormat='mm-yyyy';
            %datetick('x','mm-yyyy');
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