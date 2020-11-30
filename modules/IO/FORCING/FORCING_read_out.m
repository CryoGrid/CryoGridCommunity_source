%========================================================================
% CryoGrid FORCING class FORCING_seb

% S. Westermann, November 2020
%========================================================================

classdef FORCING_read_out
    
    properties
        forcing_index
        DATA            % forcing data time series
        TEMP            % forcing data interpolated to a timestep
        PARA            % parameters
        STATUS         
    end
    
    
    methods
        
        %constructor
        function self = FORCING_read_out(varargin)               % Temporary definition, to allow old code to run
        %function self = FORCING_seb(index, pprovider, fprovider)      % Definition to be used when old code is no longer supported
            % CONSTRUCTOR for FORCING_seb
            %   Reads in forcing data from the specified file.
            %
            %   ARGUMENTS:
            %   index:      user defined class index
            %   pprovider:  instance of PARAMETER_PROVIDER class
			%   fprovider:  instance of FORCING_PROVIDER class
            
            % The following is only needed to allow legacy code to run
            % May be removed when deprecated functions are removed
            if nargin==3
                index = varargin{1};
                pprovider = varargin{2};
				fprovider = varargin{3};
            else
                st = dbstack;
                warning(['DEPRECATION WARNING: Instantiating ' st.name '() with no arguments is deprecated.' newline,...
                         'You should update your code to take advantage of new IO interface.']);
                return
            end
            % End allow legacy code
            
            self.forcing_index = index;
            self = self.initialize();
            self = self.populate_PARA(pprovider);
            self = self.populate_DATA(fprovider);
            self = self.finalize_setup();
        end

        
        function self = initialize(self)
            % INITIALIZE  Initializes all properties needed by the class.

            self.STATUS = 0;
            self = self.initialize_PARA();
            self = self.initialize_TEMP();
            self = self.initialize_DATA();
        end
            
        
        function self = initialize_PARA(self)         
            % INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.  
            self.PARA.filename = [];
            self.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            self.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            
            self.PARA.latitude = [];  % latitude
            self.PARA.longitude = []; % longitude
            self.PARA.altitude = [];  % elevation above sea level [m]
            self.PARA.domain_depth = []; % total depth of the model domain [m]
            self.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            self.PARA.airT_height = [];  % height above ground at which air temperature (and wind speed!) from the forcing data are applied.
            self.PARA.area = [];         % area of the grid cell [m2] 
        end
        
        
        function self = initialize_TEMP(self)
            % INITIALIZE_TEMP  Initializes TEMP structure.
            
        end
        
        
        function self = initialize_DATA(self)
            % INITIALIZE_DATA  Initializes DATA structure.
            
        end
        
        
        function self = populate_PARA(self, pprovider)
            % POPULATE_PARA  Updates the PARA structure with values from pprovider.
            %
            %   ARGUMENTS:
            %   pprovider:  instance of PARAMETER_PROVIDER class
            
            self.PARA = pprovider.populate_struct(self.PARA, 'FORCING', mfilename('class'), self.forcing_index);
            if ~isfield(self.PARA, 'filepath')
                % if path to forcing data is not specified in PARA
                % assume it is in a forcing sub-directory to the main
                % script
                self.PARA.filepath = './forcing';
            end
        end
        

        function self = populate_DATA(self, fprovider)
            % POPULATE_DATA  Updates the DATA structure with values from fprovider.
            %
            %   ARGUMENTS:
			%	fprovider:	instance of FORCING class

%            self.DATA = fprovider.populate_struct(self.DATA);
 %           self.DATA.timeForcing = self.DATA.t_span;
        end
            
        
        function self = finalize_setup(self)
            
            if isempty(self.PARA.start_time) || ~ischar(self.PARA.start_time)
                self.PARA.start_time = self.DATA.timeForcing(1,1);
            else
                self.PARA.start_time = datenum(self.PARA.start_time, 'dd.mm.yyyy');
            end
            if isempty(self.PARA.end_time) || ~ischar(self.PARA.end_time)
                self.PARA.end_time = floor(self.DATA.timeForcing(end,1));
            else
                self.PARA.end_time = datenum(self.PARA.end_time, 'dd.mm.yyyy');
            end
            
            
        end


        function self = interpolate_forcing(t, self)
        
        end
        
        function self = load_forcing_from_mat(self)
            
        end
        
    end
end