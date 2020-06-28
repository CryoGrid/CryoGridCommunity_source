% base class build a model tile

classdef TILE_BUILDER
    properties
        pprovider
        cprovider
        fprovider
        
        forcing_id = 1
        grid_id = 1
        out_id = 1
        strat_linear_id =1
        strat_layers_id = 1
        strat_classes_id = 1
        %stratigraphy_id = [1 2 3]
        
        % CryoGrid class instances
        forcing
        grid
        out
        strat_linear
        strat_layers
        strat_classes
        TOP % Indicates current top position
        BOTTOM % Indicates current bottom position
        TOP_CLASS % Uppermost ground/snow class in the soil column
        BOTTOM_CLASS % Lowermost ground class in the soil column
        
        % Secondary variables: stratigraphy assembly
        stratigraphy_list
        class_list
        % etc
    end
    
    
    methods
        
        function self = TILE_BUILDER(pprovider, cprovider, fprovider, varargin)% To add to the argin: fprovider,
            % CONSTRUCTOR for TILE_BUILDER
            %   Builds a model tile, by instantiating forcing, grid and
            %   out-classes, and building the stratigraphy.
            %
            %   ARGUMENTS:
            %   pprovider:  instance of PARAMETER_PROVIDER class
            %   cprovider:  instance of CONSTANT_PROVIDER class
            %   fprovider:  instance of FORCING_PROVIDER class
            %
            %   OPTIONAL ARGUMENTS:
            %   forcing_id:       id of the forcing definition to use
            %   grid_id:          id of the grid definition to use
            %   out_id:           id of the out definition to use
            %   strat_linear_id:  id of the strat_linear definition to use
            %   strat_layers_id:  id of the strat_layers definition to use
            %   strat_classes_id: id of the strat_classes definition to use
            %
            %   User must provide the optional arguments as name,value
            %   pairs, uisng the call signature:
            %
            %   tile = TILE_BUILDER(pprovider, cprovider, fprovider, ...
            %                                  'forcing_id', 1, ...
            %                                  'grid_id', 2, ...
            %                                  'out_id', 1, ...
            %                                  'strat_linear_id', 1, ...
            %                                  'strat_layers_id',1, ...
            %                                  'strat_classes_id',1);
            %
            %   None, all or any of the id's can be specified.
            %
            %   NOTICE: id 1 refers to the first class registered in the
            %     relevant section, it does not refer to the class index!
            %     (since you may have multiple classes with index 1 listed in
            %     each section, provided they are different types of classes.)
            %
            %   We may want to be able to select by index also... currently
            %   not implemented.
            
            self.pprovider = pprovider;    % Store reference to property_provider instance
            self.cprovider = cprovider;
            self.fprovider = fprovider;
            
            % Support name-value pair arguments when constructing object
            % if the argument name is in the properties of the class
            % set the appropriate property to the argument value
            
            % Necessary to adapt the loop to the number of input arguments
            for i=1:2:nargin-3
                if any(strcmp(properties(self), varargin{i}))
                    self.(varargin{i}) = varargin{i+1};
                end
            end
            
            self = self.build_forcing();
            self = self.build_grid();
            self = self.build_out();
            self = self.build_strat_linear();
            self = self.build_strat_layers();
            self = self.build_strat_classes();
            self = self.build_stratigraphy();
            self = self.assemble_stratigraphy();
            self = set_top_depth_rel2surface(self); %function added Sebastian, required for regridding
            self = self.assemble_interactions();
            if ~isnan(self.strat_classes.snow_class.classname)
                self = self.add_SNOW(); %added Sebastian, replaces add_CHILD_snow()
            end
            if ~isempty(self.strat_classes.sleeping_classes)
                self = self.add_sleeping_classes(); %added Sebastian, appends all classes defined in STRAT_CLASSES to TOP
            end
            %          self = self.add_CHILD_snow();
            
            % continue here to implement the other classes that needs to be
            % build/instantiated...
            
            % ... others?...
            
        end
        
        
        function self = build_forcing(self)
            % BUILD_FORCING  Instantiates forcing class.
            
            [name, index] = self.pprovider.get_class_name_and_index_by_id('FORCING', self.forcing_id);
            
            class_handle = str2func(name);
            self.forcing = class_handle(index, self.pprovider, self.fprovider);
        end
        
        
        function self = build_grid(self)
            % BUILD_GRID  Instantiates grid class.
            
            [name, index] = self.pprovider.get_class_name_and_index_by_id('GRID', self.grid_id);
            
            class_handle = str2func(name);
            self.grid = class_handle(index, self.pprovider, self.forcing);
        end
        
        
        function self = build_out(self)
            % BUILD_OUT  Instantiates out class.
            
            [name, index] = self.pprovider.get_class_name_and_index_by_id('OUT', self.out_id);
            
            class_handle = str2func(name);
            self.out = class_handle(index, self.pprovider, self.forcing);
        end
        
        function self = build_strat_linear(self)
            % BUILD_STRAT_LINEAR Instantiates strat_linear class by setting the initial conditions and by linearly interpolating between the provided values.
            
            [name, index] = self.pprovider.get_class_name_and_index_by_id('STRAT_linear', self.strat_linear_id);
            
            class_handle = str2func(name);
            self.strat_linear = class_handle(index, self.pprovider, self.grid);
        end
        
        function self = build_strat_layers(self)
            % BUILD_STRAT_LAYERS  Instantiates strat_layers class from ground layers' constant properties.
            
            [name, index] = self.pprovider.get_class_name_and_index_by_id('STRAT_layers', self.strat_layers_id);
            
            class_handle = str2func(name);
            self.strat_layers = class_handle(index, self.pprovider, self.grid);
        end
        
        function self = build_strat_classes(self)
            % BUILD_STRAT_CLASSES  Instantiates strat_classes class with ground and snow
            %   modules that will be used by the model.
            
            [name, index] = self.pprovider.get_class_name_and_index_by_id('STRAT_classes', self.strat_classes_id);
            
            class_handle = str2func(name);
            self.strat_classes = class_handle(index, self.pprovider, self.grid);
        end
        
        function self = build_stratigraphy(self)
            % BUILD_STRATIGRAPHY  Builds and assembles stratigraphy from initialized ground and snow
            %   classes.
            
            % Necessary to adapt the following lines and previous functions to the case where not all the strat classes are used in a model run.
            % (Only strat_classes is supposed to be mandatory. How does the code behaves when only strat_classes is defined in the configuration file ?)
            
            self.stratigraphy_list = {self.strat_layers self.strat_layers.strat_layers_index; ...
                self.strat_classes self.strat_classes.strat_classes_index; ...
                self.strat_linear self.strat_linear.strat_linear_index};
            
            % Only initializes the classes listed in STRAT_classes class.
            % (In the legacy code, all the classes listed in the CLASS
            % section are stored and initialized in class_list)
            
            for i = 1:length(self.strat_classes.class_name)
                name = self.strat_classes.class_name{i};
                index = self.strat_classes.class_index(i);
                class_handle = str2func(name);
                self.class_list{i,1} = class_handle(index, self.pprovider, self.cprovider, self.forcing);
                self.class_list{i,2} = index;
            end
            %MODIFIED SEBASTIAN: snow class, now optional
            name = self.strat_classes.snow_class.classname;
            if ~isnan(name)
                index = self.strat_classes.snow_class.index;
                class_handle = str2func(name);
                self.class_list{size(self.class_list,1)+1,1} = class_handle(index, self.pprovider, self.cprovider, self.forcing);
                self.class_list{size(self.class_list,1),2} = index;
            end
            %MODIFIED SEBASTIAN: add sleeping classes to class_list
            if ~isempty(self.strat_classes.sleeping_classes)
                for i=1:size(self.strat_classes.sleeping_classes.class_name,1)
                    name = self.strat_classes.sleeping_classes.class_name{i,1};
                    index = self.strat_classes.sleeping_classes.class_index{i,1};
                    class_handle = str2func(name);
                    self.class_list{size(self.class_list,1)+1,1} = class_handle(index, self.pprovider, self.cprovider, self.forcing);
                    self.class_list{size(self.class_list,1),2} = index;
                end
            end

        end
        
        function self = assemble_stratigraphy(self)
            % ASSEMBLE_STRATIGRAPHY  Finds STRAT_classes with index 1 and appends  information from all other
            %   classes with index 1 to a variable_list. Indicates current top and bottom positions, finalizes STATVAR initialization and connects the stacked ground (and snow) modules with pointers (NEXT and PREVIOUS).
            
            
            for i=1:size(self.stratigraphy_list,1)  %find STRAT_class in the list
                if self.stratigraphy_list{i,2}==1
                    if strcmp(class(self.stratigraphy_list{i,1}), 'STRAT_classes')
                        class_stratigraphy = self.stratigraphy_list{i,1};
                    else
                        self.grid.variable_names = [self.grid.variable_names self.stratigraphy_list{i,1}.variable_names];
                        self.grid.variable_gridded = [self.grid.variable_gridded self.stratigraphy_list{i,1}.variable_gridded];
                    end
                end
            end
            
            
            i=1;
            for j=1:size(self.class_list,1)
                if strcmp(class(self.class_list{j,1}), class_stratigraphy.class_name{i,1}) && self.class_list{j,2}==class_stratigraphy.class_index(i,1)
                    self.TOP_CLASS = copy(self.class_list{j,1}); %make an identical copy of the class stored in class_list -> classes in class_list are fully independet of the ones in class_list
                    self.TOP_CLASS = initialize_STATVAR_from_file(self, self.TOP_CLASS, self.grid, self.forcing, class_stratigraphy.depth(i,:));  %CHANGED SEBASTIAN
                    self.TOP_CLASS = finalize_init(self.TOP_CLASS, self.forcing);   %CHANGED SEBASTIAN
                    CURRENT = self.TOP_CLASS;
                end
            end
            
            for i=2:size(class_stratigraphy.class_name,1)  %WRONG???
                for j=1:size(self.class_list,1)
                    if strcmp(class(self.class_list{j,1}), class_stratigraphy.class_name{i,1}) && self.class_list{j,2}==class_stratigraphy.class_index(i,1)
                        CURRENT.NEXT = copy(self.class_list{j,1}); %make an identical copy of the class stored in class_list
                        CURRENT.NEXT = initialize_STATVAR_from_file(self, CURRENT.NEXT, self.grid, self.forcing, class_stratigraphy.depth(i,:));  %CHANGED SEBASTIAN
                        CURRENT.NEXT = finalize_init(CURRENT.NEXT, self.forcing);  %CHANGED SEBASTIAN
                        CURRENT.NEXT.PREVIOUS = CURRENT;
                        CURRENT = CURRENT.NEXT;
                    end
                end
            end
            
            self.BOTTOM_CLASS = CURRENT;
            self.BOTTOM=Bottom(); % Instance of the Bottom class (indicates bottom position)
            self.BOTTOM = init_bottom(self.BOTTOM, self.BOTTOM_CLASS);
            self.BOTTOM_CLASS.NEXT = self.BOTTOM;
            self.TOP=Top(); % Instance of the Top class (indicates top position)
            self.TOP = init_top(self.TOP, self.TOP_CLASS);
            self.TOP_CLASS.PREVIOUS = self.TOP;
        end
        
        function self = set_top_depth_rel2surface(self)  %function required for regridding
            CURRENT = self.TOP_CLASS;
            CURRENT.STATVAR.top_depth_rel2groundSurface = 0; %set initial surface to zero
            
            CURRENT.PARA.target_grid = self.grid.GRID;
            CURRENT.PARA.target_layerThick = self.grid.LAYERTHICK;
            while ~isequal(CURRENT.NEXT, self.BOTTOM_CLASS.NEXT)
                CURRENT.NEXT.STATVAR.top_depth_rel2groundSurface = CURRENT.STATVAR.top_depth_rel2groundSurface + sum(CURRENT.STATVAR.layerThick,1);
                
                CURRENT.NEXT.PARA.target_grid = self.grid.GRID;
                CURRENT.NEXT.PARA.target_layerThick = self.grid.LAYERTHICK;
                
                CURRENT = CURRENT.NEXT;
            end
        end
        
        function self = assemble_interactions(self)
            % ASSEMBLE_INTERACTIONS  Defines interactions between model classes. Connections are established between ground modules and interactions using the pointers NEXT and PREVIOUS.
            
            CURRENT = self.TOP_CLASS;
            
            while ~isequal(CURRENT.NEXT, self.BOTTOM_CLASS.NEXT)
                ia_class = get_IA_class(class(CURRENT), class(CURRENT.NEXT));
                CURRENT.IA_NEXT = ia_class;
                CURRENT.IA_NEXT.PREVIOUS = CURRENT;
                CURRENT.IA_NEXT.NEXT = CURRENT.NEXT;
                CURRENT.NEXT.IA_PREVIOUS = CURRENT.IA_NEXT;
                CURRENT = CURRENT.NEXT;
            end
        end
        
        function self = add_SNOW(self)  %replaces add_CHILD_snow, stores an initialized SNOW object in Top(), from where it can be accessed by the stratigraphy when needed
            for i=1:size(self.stratigraphy_list,1)  %find STRAT_class in the list
                if self.stratigraphy_list{i,2}==1
                    if strcmp(class(self.stratigraphy_list{i,1}), 'STRAT_classes')
                        class_stratigraphy = self.stratigraphy_list{i,1};
                    end
                end
            end
            
            for i=1:size(self.class_list,1)  %find snow in the class list
                if strcmp(class(self.class_list{i,1}), class_stratigraphy.snow_class.classname) && self.class_list{i,2} == class_stratigraphy.snow_class.index
                    self.TOP.STORE.SNOW=copy(self.class_list{i,1});
                end
            end
            self.TOP.STORE.SNOW = finalize_init(self.TOP.STORE.SNOW, self.forcing);
            
        end
        
        function self = add_sleeping_classes(self)  %adds sleeping_classes to TOP.STORE.SLEEPING
            for i=1:size(self.stratigraphy_list,1)  %find STRAT_class in the list
                if self.stratigraphy_list{i,2}==1
                    if strcmp(class(self.stratigraphy_list{i,1}), 'STRAT_classes')
                        class_stratigraphy = self.stratigraphy_list{i,1};
                    end
                end
            end
            
            for i=1:size(self.class_list,1)  %find sleeping_classes in the class list
                for j = 1:size(class_stratigraphy.sleeping_classes.class_name,1)
                    if strcmp(class(self.class_list{i,1}), class_stratigraphy.sleeping_classes.class_name{j,1}) && self.class_list{i,2} == class_stratigraphy.sleeping_classes.class_index{j,1}
                        self.TOP.STORE.SLEEPING{j,1} = copy(self.class_list{i,1});
                        self.TOP.STORE.SLEEPING{j,1} = finalize_init(self.TOP.STORE.SLEEPING{j,1}, self.forcing);
                        self.TOP.STORE.SLEEPING{j,2} = class_stratigraphy.sleeping_classes.class_index{j,1};
                    end
                end
            end            
        end

        %         function self = add_CHILD_snow(self)  %discontinued
        %             % ADD_CHILD_SNOW  Connects the uppermost ground module to a snow module that is built as a child, using an interaction. The snow child class is at first initilized with empty pointers and zero snow conditions.
        %
        %
        %
        %             %replace by matrix
        %
        %             if strcmp(class(self.TOP_CLASS), 'GROUND_freeW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_bucketW')
        %                 self.TOP_CLASS.IA_CHILD = IA_SNOW_GROUND();
        %             elseif strcmp(class(self.TOP_CLASS), 'GROUND_freeW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_crocus')
        %                 self.TOP_CLASS.IA_CHILD = IA_SNOW_GROUND_crocus();
        %             elseif strcmp(class(self.TOP_CLASS), 'GROUND_fcSimple_salt_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_bucketW')
        %                 self.TOP_CLASS.IA_CHILD = IA_SNOW_GROUND_fcSimple_salt();
        %             elseif strcmp(class(self.TOP_CLASS), 'GROUND_fcSimple_salt_seb_snow') && strcmp(class(snow), 'SNOW_crocus_no_inheritance')
        %                 self.TOP_CLASS.IA_CHILD = IA_SNOW_GROUND_fcSimple_salt_crocus();
        %             elseif strcmp(class(self.TOP_CLASS), 'GROUND_freeW_bucketW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_bucketW')
        %                 self.TOP_CLASS.IA_CHILD = IA_SNOW_GROUND();
        %             elseif strcmp(class(self.TOP_CLASS), 'GROUND_freeW_bucketW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_crocus')
        %                 self.TOP_CLASS.IA_CHILD = IA_SNOW_GROUND_crocus();
        %             elseif strcmp(class(self.TOP_CLASS), 'GROUND_freeW_bucketW_seb_snow') && strcmp(class(snow), 'SNOW_crocus_no_inheritance')
        %                 self.TOP_CLASS.IA_CHILD = IA_SNOW_GROUND_crocus();
        %             elseif strcmp(class(self.TOP_CLASS), 'GROUND_freeW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb')
        %                 self.TOP_CLASS.IA_CHILD = IA_SNOW_GROUND();
        %             end
        %
        %
        %             CURRENT = self.TOP_CLASS.IA_CHILD; %change to interaction class
        %             CURRENT.STATUS = 0; %snow initially inactive
        %             CURRENT.IA_PARENT_GROUND = self.TOP_CLASS;
        %             CURRENT.IA_CHILD_SNOW = snow;
        %             CURRENT.IA_CHILD_SNOW = initialize_zero_snow(CURRENT.IA_CHILD_SNOW, CURRENT.IA_PARENT_GROUND);
        %
        %             self.TOP_CLASS.IA_CHILD = CURRENT; %reassign to ground
        %         end
        
        
        function ground = initialize_STATVAR_from_file(self, ground, grid, forcing, depths);
            variables = fieldnames(ground.STATVAR);
            range = (grid.MIDPOINTS > depths(1,1) & grid.MIDPOINTS <= depths(1,2));
            ground.STATVAR.layerThick = grid.LAYERTHICK(range,1);
            ground.STATVAR.upperPos = forcing.PARA.altitude - depths(1,1);
            ground.STATVAR.lowerPos = forcing.PARA.altitude - depths(1,2);
            
            for j=1:size(variables,1)
                size(grid.variable_names,2);
                for i=1:size(grid.variable_names,2)
                    if strcmp(variables{j,1}, grid.variable_names{1,i})
                        variables{j,1};
                        range;
                        
                        ground.STATVAR.(variables{j,1}) = grid.variable_gridded(range,i);
                    end
                end
            end
            
        end
        
        %        function self = build_stratigraphy(self)
        %            % BUILD_STRATIGRAPHY 1st method used to build the stratigraphy and instantiate the ground and snow classes
        %            % To remove if the method currently implemented is functional !
        
        %            [name, index] = self.pprovider.get_class_name_and_index_by_id('STRATIGRAPHY', self.stratigraphy_id);
        
        %             if name == 'STRAT_classes'
        %               for j = 1:length(self.stratigraphy{i,1}.class_name)
        %                   name = self.stratigraphy{i,1}.class_name{j};
        %                   index = self.stratigraphy{i,1}.class_index(j);
        %                   class_handle = str2func(name);
        %                   self.class_list{j,1} = class_handle(index, self.pprovider, self.grid)
        %                   self.class_list{j,2} = index;
        %               end
        %               name = self.stratigraphy{i,1}.snow_class.name;
        %               index = self.stratigraphy{i,1}.snow_class.index;
        %               class_handle = str2func(name);
        %               self.class_list{j+1,1} = class_handle(index, self.pprovider, self.grid)
        %               self.class_list{j+1,2} = index;
        %             end
        %             error([mfilename('class') '.build_stratigraphy() method is not yet implemented.']);
        %             for i = 1:length(self.stratigraphy{2,1}.class_name)
        %               get_class_id_by_name_and_index(self, section, name, index)
        % %           end
        %         end
    end
end
