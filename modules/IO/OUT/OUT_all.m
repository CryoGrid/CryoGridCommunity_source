% forcing data


classdef OUT_all

    properties
		out_index
        STRATIGRAPHY
        TIMESTAMP
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
		
	end
    
    
    methods
		
		function self = OUT_all(varargin)               % Temporary definition, to allow old code to run
        %function self = OUT_all(index, pprovider, forcing)      % Definition to be used when old code is no longer supported
            % CONSTRUCTOR for OUT_all
            %   Reads in out data from the specified file.
            %
            %   ARGUMENTS:
            %   index:      user defined class index
            %   pprovider:  instance of PARAMETER_PROVIDER class
			%	forcing:	instance of FORCING class
            
            % The following is only needed to allow legacy code to run
            % May be removed when deprecated functions are removed
            if nargin==3
                index = varargin{1};
                pprovider = varargin{2};
				forcing = varargin{3};
            else
                st = dbstack;
                warning(['DEPRECATION WARNING: Instantiating ' st.name '() with no arguments is deprecated.' newline,...
                         'You should update your code to take advantage of new IO interface.']);
                return
            end
            % End allow legacy code
            
            self.out_index = index;
            self = self.initialize();
            self = self.populate_PARA(pprovider);
            self = self.finalize_setup(forcing);
        end
		
		function self = initialize(self)
            % INITIALIZE  Initializes all properties needed by the class.
			
			self.STRATIGRAPHY = [];
			self.TIMESTAMP = [];
			self.MISC = [];
			self.OUTPUT_TIME = [];
			self.SAVE_TIME = [];
            self = self.initialize_PARA();
            self = self.initialize_TEMP();
        end
        
        function self = initialize_PARA(self)         
            % INITIALIZE_PARA  Initializes PARA structure.

            self.PARA.output_timestep = [];
            self.PARA.save_date = [];
            self.PARA.save_interval = [];
        end
		
	    function self = initialize_TEMP(self)
            % INITIALIZE_TEMP  Initializes TEMP structure.
            
            self.TEMP = struct();
        end	
		
		function self = populate_PARA(self, pprovider)
            % POPULATE_PARA  Updates the PARA structure with values from pprovider.
            %
            %   ARGUMENTS:
            %   pprovider:  instance of PARAMETER_PROVIDER class
            
            self.PARA = pprovider.populate_struct(self.PARA, 'OUT', mfilename('class'), self.out_index);
        end
		
		function self = finalize_setup(self, forcing)
			% FINALIZE_SETUP  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
			
			%	ARGUMENTS:
			%	forcing:	instance of FORCING class
			
			self.OUTPUT_TIME = forcing.PARA.start_time + self.PARA.output_timestep;
            if isempty(self.PARA.save_interval) || isnan(self.PARA.save_interval) 
                self.SAVE_TIME = forcing.PARA.end_time;
            else
                self.SAVE_TIME = min(forcing.PARA.end_time,  datenum([self.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + self.PARA.save_interval)], 'dd.mm.yyyy'));
            end
		end
		
		function self = store_OUT(self, t, TOP_CLASS, BOTTOM, forcing, run_number, timestep, result_path)
            
            if t==self.OUTPUT_TIME
                %if id == 1
                disp([datestr(t)])
                %end
                %labBarrier
                self.TIMESTAMP=[self.TIMESTAMP t];
                
                %self.STRATIGRAPHY{1,size(self.STRATIGRAPHY,2)+1} = copy(TOP_CLASS);  %append new stratigraphy, should be made more sophisticated by not adding instaneous values, but averaging/accumulating variables
                %self.STRATIGRAPHY{1,size(self.STRATIGRAPHY,2)+1} = [TOP_CLASS.STATVAR.T; TOP_CLASS.NEXT.STATVAR.T];  
                CURRENT =TOP_CLASS;
                if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
                    self.MISC=[self.MISC [CURRENT.CHILD.STATVAR.T(1,1); CURRENT.CHILD.STATVAR.layerThick(1,1)]]; 
                else
                    self.MISC=[self.MISC [NaN; NaN]];
                end
                result={};
                while ~isequal(CURRENT, BOTTOM)
                    if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
                        res=copy(CURRENT.CHILD);
                        res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_NEXT=[];  res.PARENT = []; %cut all dependencies
                        result=[result; {res}];
                    end
                    res = copy(CURRENT);
                    if isprop(res, 'LUT')
                        res.LUT =[];  %remove look-up tables, runs out of memeory otherwise
                    end
                    res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_NEXT=[];  %cut all dependencies
                    if isprop(res, 'CHILD')
                        res.CHILD = [];
                        res.IA_CHILD =[];
                    end
                    result=[result; {res}];
                    CURRENT = CURRENT.NEXT;
                end
                self.STRATIGRAPHY{1,size(self.STRATIGRAPHY,2)+1} = result;
                
                self.OUTPUT_TIME = self.OUTPUT_TIME + self.PARA.output_timestep;
                if t==self.SAVE_TIME 
                   if ~(exist([result_path run_number])==7)
                       mkdir([result_path run_number])
                   end
                   save([result_path run_number '/' run_number '_' datestr(t,'yyyy') '.mat'], 'self')
                   self.STRATIGRAPHY=[];
                   self.TIMESTAMP=[];
                   self.MISC=[];
                   self.SAVE_TIME = min(forcing.PARA.end_time,  datenum([self.PARA.save_date num2str(str2num(datestr(self.SAVE_TIME,'yyyy')) + self.PARA.save_interval)], 'dd.mm.yyyy'));
                end
            end
        end

        function xls_out = write_excel(self)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			
            xls_out = {'OUT','index',NaN,NaN;'OUT_all',1,NaN,NaN;'output_timestep',0.250000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'OUT_END',NaN,NaN,NaN};
        end
         
        
		% ==========================================
        % DEPRECATED METHODS
        % to be deleted when new implementation
        % is validated for backwards compatibility
        % ==========================================
		
		% function res = initialize_out_all(run_info) %discontinued
%             res.PARA.output_timestep = 1/4;
%             res.PARA.save_date = '01.09.';
%             res.PARA.save_interval = 1;
%             res.OUTPUT_TIME = run_info.START_TIME + res.PARA.output_timestep;
%             if isempty (res.PARA.save_interval)
%                 res.SAVE_TIME = run_info.END_TIME;
%             else
%                 res.SAVE_TIME = min(run_info.END_TIME,  datenum([res.PARA.save_date num2str(str2num(datestr(run_info.START_TIME,'yyyy')) + res.PARA.save_interval)], 'dd.mm.yyyy'));
%             end
%         end
		
        function self = provide_variables(self)
			st = dbstack;
            warning(['DEPRECATION WARNING: Method ' st.name '() is deprecated and will be removed.' newline,...
                     'Use PARAMETER_PROVIDER class to obtain parameter values.']);
            self.PARA.output_timestep = [];
            self.PARA.save_date = [];
            self.PARA.save_interval = [];
        end
        
        function self = initalize_from_file(self, section)
			st = dbstack;
			warning(['DEPRECATION WARNING: Method ' st.name '() is deprecated and will be removed.' newline,...
                     'Use PARAMETER_PROVIDER class to obtain parameter values.']);
					 
			variables = fieldnames(self.PARA);
			for i=1:size(variables,1)
				for j=1:size(section,1)
					if strcmp(variables{i,1}, section{j,1})
						self.PARA.(variables{i,1}) = section{j,2};
					end
				end
			end
		end

		function self = initialize_from_ini(self, ini)
			% INITIALIZE_FROM_INI  Initializes the variables from output structure of the ini parser, and compares the
			%	names of the variables from the class to the ini structure.
			% 	If the variables from the class mismatch the variables from
			% 	the ini file, an error message is displayed.
			
			%	ARGUMENTS:
			%	ini:	output structure from the ini parser
					
			st = dbstack;
            warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
                     'Code should be moved to new PARAMETER_PROVIDER class ',...
                     'to streamline file access and the population of parameters.']);
					 
			ini_variables = fields(ini.OUT_all);
			out_variables = fieldnames(self.PARA);
			ismatch_class_ini_variables(ini_variables, out_variables) 
			for i=1:length(out_variables)
				for j=1:length(ini_variables)
					if strcmp(out_variables{i,1},ini_variables{i,1})
						self.PARA.(out_variables{i,1}) = ini.out_all.(ini_variables{i,1})
					end
				end
			end
		end
			
		% 
        function self = complete_init_out(self, forcing)
			% COMPLETE_INIT_OUT  Completes the initialization of the member variables of the class by creating them based on out.PARA and forcing.PARA.
			
			%	ARGUMENTS:
			%	forcing:	instance of FORCING class
			
			st = dbstack;
            warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
            'Parameter initialization should be finalized in the ' mfilename('class') '.finalize_setup() ']);
            
            self.OUTPUT_TIME = forcing.PARA.start_time + self.PARA.output_timestep;
            if isempty(self.PARA.save_interval) || isnan(self.PARA.save_interval) 
                self.SAVE_TIME = forcing.PARA.end_time;
            else
                self.SAVE_TIME = min(forcing.PARA.end_time,  datenum([self.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + self.PARA.save_interval)], 'dd.mm.yyyy'));
            end
        end 
        
    end
end