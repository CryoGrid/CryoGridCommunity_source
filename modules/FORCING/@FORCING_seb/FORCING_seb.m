%========================================================================
% CryoGrid FORCING class FORCING_seb
% simple model forcing for GROUND classes computing the surface energy balance 
% (keyword “seb”). The data must be stored in a Matlab “.mat” file which contains 
% a struct FORCING with field “data”, which contain the time series of the actual 
% forcing data, e.g. FORCING.data.Tair contains the time series of air temperatures. 
% Have a look at the existing forcing files in the folder “forcing” and prepare 
% new forcing files in the same way. The mandatory forcing variables are air temperature 
% (Tair, in degree Celsius), incoming long-wave radiation (Lin, in W/m2), 
% incoming short-.wave radiation (Sin, in W/m2), absolute humidity (q, in 
% kg water vapor / kg air), wind speed (wind, in m/sec), rainfall (rainfall, in mm/day), 
% snowfall (snowfall, in mm/day) and timestamp (t_span, 
% in Matlab time / increment 1 corresponds to one day). 
% IMPORTANT POINT: the time series must be equally spaced in time, and this must be 
% really exact. When reading the timestamps from an existing data set (e.g. an Excel file),
% rounding errors can result in small differences in the forcing timestep, often less 
% than a second off. In this case, it is better to manually compile a new, equally spaced 
% timestep in Matlab.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef FORCING_seb
    
    properties
        forcing_index
        DATA            % forcing data time series
        TEMP            % forcing data interpolated to a timestep
        PARA            % parameters
        STATUS         
    end
    
    
    methods
        
        %constructor
        function self = FORCING_seb(varargin)               % Temporary definition, to allow old code to run
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

            self.PARA.filename = [];   %filename of Matlab file containing forcing data
            self.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            self.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            self.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            self.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
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
            
            self.TEMP.snowfall=0;
            self.TEMP.rainfall=0;
            self.TEMP.Lin=0;
            self.TEMP.Sin=0;
            self.TEMP.Tair=0;
            self.TEMP.wind=0;
            self.TEMP.q=0;
            self.TEMP.p=0;
			%self.TEMP.t=[];
        end
        
        
        function self = initialize_DATA(self)
            % INITIALIZE_DATA  Initializes DATA structure.
            
            self.DATA.rainfall    = [];
            self.DATA.snowfall    = [];
            self.DATA.Tair        = [];
            self.DATA.Lin         = [];
            self.DATA.Sin         = [];
            self.DATA.q           = [];
            self.DATA.wind        = [];
            self.DATA.timeForcing = [];  % timeForcing = t_span
            self.DATA.t_span      = [];
            self.DATA.z           = [];
            self.DATA.lon         = [];
            self.DATA.lat         = [];
            self.DATA.p           = [];
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

            self.DATA = fprovider.populate_struct(self.DATA);
            self.DATA.timeForcing = self.DATA.t_span;
        end
            
        
        function self = finalize_setup(self)
            % FINALIZE_SETUP  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
            
            % apply rain and snow fractions
            self.DATA.rainfall = self.DATA.rainfall.*self.PARA.rain_fraction;
            self.DATA.snowfall = self.DATA.snowfall.*self.PARA.snow_fraction;        
            
            % Update spatial data if included in forcing DATA
            if ~isempty(self.DATA.z)
                self.PARA.altitude = round(slef.DATA.z);
            end
            if ~isempty(self.DATA.lon) && ~isempty(self.DATA.lat)
                 self.PARA.longitude = self.DATA.lon;
                 self.PARA.latitude  = self.DATA.lat;
            end

            if std(self.DATA.timeForcing(2:end,1)-self.DATA.timeForcing(1:end-1,1))>=1e-10
                error('timestamp of forcing data is not in regular intervals -> check, fix and restart')
                %self.STATUS=0;
                %return
            %else
            %    self.STATUS=1;
            end
            
            self.STATUS = 1;

            % here, consistency checks, RH->q calculation, set threhsolds for wind, etc. could be placed

            % forcing.DATA.wind(forcing.DATA.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence
            self.DATA.Lin(find(self.DATA.Lin==0)) = 5.67e-8 .* (self.DATA.Tair(find(self.DATA.Lin==0))+273.15).^4;

            %set pressure to mean pressure at corresponding altitude (international
            %altitude formula) if now provided by the forcing time series
            if isempty(self.DATA.p)
                self.DATA.p = self.DATA.Tair.*0 + 1013.25.*100.*(1-0.0065./288.15.*self.PARA.altitude).^5.255;
            end

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

        function xls_out = write_excel(self)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			
            xls_out = {'FORCING','index',NaN,NaN;'FORCING_seb',1,NaN,NaN;NaN,NaN,NaN,NaN;'filename',NaN,NaN,NaN;'start_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the first timestamp of the forcing data set will be used';'end_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the last timestamp of the forcing data set will be used';'rain_fraction',1,'[-]','rainfall in forcing file multiplied by this number';'snow_fraction',1,'[-]','snowfall in forcing file multiplied by this number';'latitude',NaN,'[degree]','geographical coordinates';'longitude',NaN,'[degree]',NaN;'altitude',NaN,'[m]','a.s.l.';'domain_depth',100,'[m]','should match a GRID point, model domain extends to this depth';'heatFlux_lb',0.0500000000000000,'[W/m2]','geothermal heat flux';'airT_height',2,'[m]','height of air temperature';'FORCING_END',NaN,NaN,NaN};
        end
        
        
%         % ==========================================
%         % DEPRECATED METHODS
%         % to be deleted when new implementation
%         % is validated for backwards compatibility
%         % ==========================================
%         
%         
%         function self = initalize_from_file(self, section)
%             st = dbstack;
%             warning(['DEPRECATION WARNING: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Use PARAMETER_PROVIDER class to obtain parameter values.']);
%             
%             for i=1:size(section,1)
%                 if strcmp(section{i,1}, 'filename')
%                     self.PARA.filename = section{i,2};
%                 end
%                 if strcmp(section{i,1}, 'rain_fraction')
%                     self.PARA.rain_fraction = section{i,2};
%                 end
%                 if strcmp(section{i,1}, 'snow_fraction')
%                     self.PARA.snow_fraction = section{i,2};
%                 end
%                 if strcmp(section{i,1}, 'latitude')
%                     self.PARA.latitude = section{i,2};
%                 end
%                 if strcmp(section{i,1}, 'longitude')
%                     self.PARA.longitude = section{i,2};
%                 end
%                 if strcmp(section{i,1}, 'altitude')
%                     self.PARA.altitude = section{i,2};
%                 end
%                 if strcmp(section{i,1}, 'domain_depth')
%                     self.PARA.domain_depth = section{i,2};
%                 end
%                 if strcmp(section{i,1}, 'heatFlux_lb')
%                     self.PARA.heatFlux_lb = section{i,2};
%                 end
%                 if strcmp(section{i,1}, 'airT_height')
%                     self.PARA.airT_height = section{i,2};
%                 end
%                 if strcmp(section{i,1}, 'start_time')
%                     self.PARA.start_time = section{i,2};
%                 end
%                 if strcmp(section{i,1}, 'end_time')
%                     self.PARA.end_time = section{i,2};
%                 end 
%             end
%         end
%         
%         
% 		% 
%         function self = initialize_from_ini(self, ini)
% 			% INITIALIZE_FROM_INI  Initializes the variables from output structure of the ini parser, and compares the
% 			%	names of the variables from the class to the ini structure.
% 			% 	If the variables from the class mismatch the variables from
% 			% 	the ini file, an error message is displayed.
% 			
% 			%	ARGUMENTS:
% 			%	ini:	output structure from the ini parser
% 			
%             st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Code should be moved to new PARAMETER_PROVIDER class ',...
%                      'to streamline file access and the population of parameters.']);
%             
% 			ini_variables = fields(ini.FORCING_seb);
%             forcing_variables = fields(self.PARA);
%             ismatch_class_ini_variables(ini_variables, forcing_variables) 
%             for i=1:length(forcing_variables)
%                 for j=1:length(ini_variables)
%                     if strcmp(forcing_variables{i,1},ini_variables{i,1})
%                         self.PARA.(forcing_variables{i,1}) = ini.FORCING_seb.(ini_variables{i,1});
%                     end
%                 end
%             end
%         end
% 		
%         
%         function self = set_parameters(self, filename, rain_fraction, snow_fraction, altitude)
%             st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Parameters should be initialized in the ' mfilename('class') '.initialize_variables() ',... 
%                      'method and populated using the ' mfilename('class') '.populate() method.']);
%             
%             self.PARA.filename = filename;
%             self.PARA.rain_fraction = rain_fraction;
%             self.PARA.snow_fraction = snow_fraction;
%             self.PARA.altitude = altitude;
%         end
%         
    end
end