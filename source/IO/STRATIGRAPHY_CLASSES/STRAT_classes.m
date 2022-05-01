%========================================================================
% CryoGrid STRATIGRAPHY_CLASSES class STRAT_classes 
% STRAT_classes defines the initial CryoGrid stratigraphy of stratigraphy 
% classes.
% Each stratigraphy class is identified by its name and a class_index 
% (integer number)  which makes it possible to define the same stratigraphy
% classs everal times with different parameters. For stratigraphy classes 
% interacting with a snow cover, a SNOW class must be % defined.
% Sleeping classes are defined within the stratigraphy are initialized, 
% but not included in the initial CryoGrid stratigraphy. Instead, they 
% are made accessible when requested by a trigger.
% An example is the formation of a lake on initially dry ground. A LAKE
% class is initialized as sleeping class, and added on top of the
% stratigraphy when sufficient water has accumulated.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef STRAT_classes < matlab.mixin.Copyable

    
    properties
		strat_classes_index
		PARA
        CONST
        depth
        class_name
        class_index
        snow_class
        sleeping_classes
    end
    
    methods
        
% 		function self = STRAT_classes(varargin)               % Temporary definition, to allow old code to run
%         %function self = STRAT_classes(index, pprovider, grid)      % Definition to be used when old code is no longer supported
%             % CONSTRUCTOR for STRAT_classes
%             %   Reads in variable profile from the specified file.
%             %
%             %   ARGUMENTS:
%             %   index:      user defined class index
%             %   pprovider:  instance of PARAMETER_PROVIDER class
% 			%	grid:		instance of GRID class
%             
%             % The following is only needed to allow legacy code to run
%             % May be removed when deprecated functions are removed
%             if nargin==3
%                 index = varargin{1};
%                 pprovider = varargin{2};
% 				grid = varargin{3};
%             else
%                 st = dbstack;
%                 warning(['DEPRECATION WARNING: Instantiating ' st.name '() with no arguments is deprecated.' newline,...
%                          'You should update your code to take advantage of new IO interface.']);
%                 return
%             end
%             % End allow legacy code
%             
%             self.strat_classes_index = index;
%             self = self.initialize();
%             self = self.populate_variables(pprovider);
%             self = assign_sleeping_classes(self);
%             self = self.finalize_setup(grid);
%         end

%         function self = initialize_excel(self)
%             
%         end
		
% 		function self = initialize(self)
%             % INITIALIZE  Initializes all properties needed by the class.
% 
%             self.depth = [];
% 			self.class_name = [];
% 			self.class_index = [];
% 			self = self.initialize_PARA();
% 			self = self.initialize_snow_class();
%         end
		
		function self = provide_PARA(self)
			% INITIALIZE_PARA  Initializes PARA structure.
            
			self.PARA.classes = [];
			self.PARA.snow_class_name = [];
			self.PARA.snow_class_index = [];
            self.PARA.sleeping_classes_name = [];
            self.PARA.sleeping_classes_index = [];
        end
		
        function self = provide_CONST(self)

        end
        
        function self = provide_STATVAR(self)

        end       
        
        
        %-------------param file generation-----
         function stratigraphy = param_file_info(stratigraphy)
             stratigraphy = provide_PARA(stratigraphy);
             %default
             stratigraphy.PARA.default_value = [];
             
             stratigraphy.PARA.STATVAR = [];
             stratigraphy.PARA.class_category = 'STRATIGRAPHY_CLASSES';

             stratigraphy.PARA.options.classes.name = 'STRAT_MATRIX';
             stratigraphy.PARA.options.classes.entries_x = {'class_name', 'class_index'};
             stratigraphy.PARA.options.classes.entries_y = {0; 20};
             stratigraphy.PARA.comment.classes = {'stratigraphy of subsurface classes'};
             
             stratigraphy.PARA.comment.snow_class_name = {'snow class to be used in the model run'};
             
             stratigraphy.PARA.options.sleeping_classes_name.name = 'H_LIST';
             stratigraphy.PARA.options.sleeping_classes_name.entries_x = {};
             stratigraphy.PARA.comment.sleeping_classes_name = {'classes that are added to the stratigraphy during the model run'};
             
             stratigraphy.PARA.options.sleeping_classes_index.name = 'H_LIST';
             stratigraphy.PARA.options.sleeping_classes_index.entries_x = {};
         end
        
        
        
        
% 		function self = initialize_snow_class(self)
%             % INITIALIZE_SNOW_CLASS  Initializes the snow class name and index.
% 
% 			self.snow_class.classname = {};
% 			self.snow_class.index = [];
% 		end		
% 		
% 		function self = populate_variables(self, pprovider)
% 			% POPULATE_VARIABLES  Updates the PARA structure with values from pprovider. Assigns values from the PARA structure to the corresponding class properties.
%             %
%             %   ARGUMENTS:
%             %   pprovider:  instance of PARAMETER_PROVIDER class
% 			
% 			self.PARA = pprovider.populate_struct(self.PARA, 'STRAT_classes', mfilename('class'), self.strat_classes_index);
% 			
% 			fn_substruct = fieldnames(self.PARA.classes);
% 			fn = fieldnames(self.PARA);
% 			p = properties(self);
% 			for i = 1:size(fn_substruct, 1)
% 				if any(strcmp(p, fn_substruct(i)))
% 					index = find(strcmp(p, fn_substruct(i)));
% 					self.(p{index}) = self.PARA.classes.(fn_substruct{i});
% 				end
%             end
%             
%             self.depth = cell2mat(self.depth);
%             self.class_index = cell2mat(self.class_index);
%             
% 			for j = 1:size(fn, 1)
%                 if any(strncmp(fn(j),p,10))
%                     index = find(strncmp(p, fn(j),10));
% 					if isstruct(self.(p{index}))
%                         fn_struct_p = fieldnames(self.(p{index}));
%                         if endsWith(fn{j},fn_struct_p)
%                             index2 = find(~cellfun(@isempty,regexp(fn{j},fn_struct_p)));
%                             self.(p{index}).(fn_struct_p{index2}) = self.PARA.(fn{j});
%                         end
%                     end
%                 end
%             end
%         end
%         
%         function self = assign_sleeping_classes(self)  %removes the classes without depth assigned and stores the entire class list in sleeping_classes
%             pos_not_nan = ~isnan(self.depth);
%             %self.sleeping_classes = [];
% 
%             if sum(double(~pos_not_nan))>0
%                 self.sleeping_classes.class_name = self.class_name(~pos_not_nan,1);
%                 self.sleeping_classes.class_index = self.PARA.classes.class_index(~pos_not_nan,1);
%                 self.depth = self.depth(pos_not_nan,1);
%                 self.class_name = self.class_name(pos_not_nan,1);
%                 self.PARA.classes.depth = self.PARA.classes.depth(pos_not_nan,1);
%                 self.PARA.classes.class_index = self.PARA.classes.class_index(pos_not_nan,1);
%                 self.PARA.classes.class_name = self.PARA.classes.class_name(pos_not_nan,1);
%             end
%         end
% 		
% 		function self = finalize_setup(self, grid)
% 			% FINALIZE_SETUP  Performs all additional property
%             %   initializations and modifications. Checks for some (but not
%             %   all) data validity.
% 			
% 			%	ARGUMENTS:
% 			%	grid:	instance of GRID class
% 			
%             depth2=self.depth;
%             for i=1:size(depth2)-1
%                 depth2(i,1) = self.depth(i+1,1);
%             end
%             depth2(end,1) = grid.GRID(end,1);
%             self.depth = [self.depth depth2];
% 		end
% 
%         function xls_out = write_excel(self)
% 			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
% 			
%             xls_out = {'STRATIGRAPHY','index',NaN;'STRAT_classes',1,NaN;'depth - top of layer','classname','index';'[m]',NaN,NaN;'TOP',NaN,NaN;0,NaN,1;2,NaN,1;'BOTTOM',NaN,NaN;NaN,NaN,NaN;'snow_class',NaN,1;NaN,NaN,NaN;'STRATIGRAPHY_END',NaN,NaN};
%         end
        
% 		% ==========================================
%         % DEPRECATED METHODS
%         % to be deleted when new implementation
%         % is validated for backwards compatibility
%         % ==========================================
% 		
%         function self = initalize_from_file(self, section)
% 			st = dbstack;
% 			warning(['DEPRECATION WARNING: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Use PARAMETER_PROVIDER class to obtain parameter values.']);
% 					 
%             for i=1:size(section,1)
%                 if strcmp(section{i,1}, 'snow_class')
%                     self.snow_class.name = section{i,2};
%                     self.snow_class.index = section{i,3};
%                 end
%             end
%             pos_list = get_range_TOP_BOTTOM(section);
%             self.depth = cell2mat(section(pos_list(1,1):pos_list(1,2), 1));
%             self.class_name =  section(pos_list(1,1):pos_list(1,2), 2);
%             self.class_index = cell2mat(section(pos_list(1,1):pos_list(1,2), 3));
%         end
%         
% 		function self = initialize_from_table(self, table)
% 			% INITIALIZE_FROM_TABLE  Initializes the variables from the output table of the csv parser.
% 			
% 			%	ARGUMENTS:
% 			%	table:	Matlab output table from the csv parser
% 			
% 			st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Code should be moved to new PARAMETER_PROVIDER class ',...
%                      'to streamline file access and the population of parameters.']); 
% 					 
% 			self.depth = table.depth;
%             for i = 1:length(table.Properties.VariableNames)
% 				if strcmp(table.Properties.VariableNames{i}, 'snow_class')
% 					self.snow_class.name = table.snow_class{1};
%                     self.snow_class.index = table.snowindex(1);
%                 end
%             end
%             self.class_name =  table.class_name;
%             self.class_index = table.groundindex;
%         end
% 			
%         function self = interpolate_to_grid(self, grid)
% 			st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Parameter initialization should be finalized in the ' mfilename('class') '.finalize_setup() ']);
%            
%             depth2=self.depth;
%             for i=1:size(depth2)-1
%                 depth2(i,1) = self.depth(i+1,1);
%             end
%             depth2(end,1) = grid.GRID(end,1);
%             self.depth = [self.depth depth2];
%         end
        
    end
    
end

