%========================================================================
% CryoGrid FORCING class FORCING_Tair_Tlb_fixed
%
% simple model forcing for GROUND classes depending only on an upper
% boundary temperature (Tair) forcing, and a fixed lower boundary
% temperature.
%
% The data is obtained using the READ_FORCING_mat class. See this class for
% instructions about mat-file data format.
%
% The mandatory forcing variables are:
%
% Tair:      Air temperature (in degree Celsius)
% T_lb:      Lower boundary temperature (in degree Celsius)
%
% T_lb is specified in the parameter input file, not in the forcing mat 
% file.
%
% One array of timestamps must be provided (t_stamp, in Matlab time / 
% increment 1 corresponds to one day). 
%
% IMPORTANT POINT: the time series must be equally spaced in time, and this 
% must be really exact. When reading the timestamps from an existing data 
% set (e.g. an Excel file), rounding errors can result in small differences 
% in the forcing timestep, often less than a second off. In this case, it 
% is better to manually compile a new, equally spaced timestep in Matlab.
%
% T. Ingeman-Nielsen, October 2022
%========================================================================

classdef FORCING_Tair_Tlb_fixed_mat < FORCING_base & READ_FORCING_mat
    
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
 
        
        function forcing = finalize_init(forcing, tile)
            % FINALIZE_INIT  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
            
            variables = {'Tair'};
            [data, times] = read_mat([forcing.PARA.forcing_path forcing.PARA.filename], variables);
            
            for i=1:size(variables,1)
                if isfield(data, variables{i,1})
                    forcing.DATA.(variables{i,1}) = data.(variables{i,1});
                end
            end

            forcing.DATA.timeForcing = times;

            forcing = check_and_correct(forcing); % Remove known errors
            forcing = set_start_and_end_time(forcing); % assign start/end time
            
            %initialize TEMP
            forcing.TEMP.Tair = 0;
			forcing.TEMP.T_lb = forcing.PARA.T_lb;
            
        end

        function forcing = interpolate_forcing(forcing, tile)
            % Interpolate forcing data to timestep tile.t
            forcing = interpolate_forcing@FORCING_base(forcing, tile);

			forcing.TEMP.T_lb = forcing.PARA.T_lb;
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
   
end