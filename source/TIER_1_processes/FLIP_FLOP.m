%========================================================================
% CryoGrid TIER1 library class FLIP_FLOP, designed to store all STATVAR
% from a class in a model phase, and then read these STATVAR in a read phase. This makes multi-millennial runs possible for some applications. 
% S. Westermann, November 2021
%========================================================================

classdef FLIP_FLOP < BASE
    
    properties
       STORE 
       MODEL_CLASS
    end
    
    methods
        
        function ground = provide_PARA_flip_flop(ground)
            ground.PARA.number_of_model_years = [];  %the N in 1, 2, ...N
            ground.PARA.start_store_year = []; % the N_s in N_s, 6, 7, ... N
            %ground.PARA.end_store_year = [];
            ground.PARA.number_of_read_years = [];
            ground.PARA.save_date = [];
            ground.PARA.store_interval = []; %in days
            
        end
        
        function ground = provide_CONST_flip_flop(ground)
            ground.CONST.day_sec =[];
        end
        
        
        
        function ground = finalize_init_flip_flop(ground, tile)
            
                save_day = str2num(ground.PARA.save_date(1:2));
                save_month = str2num(ground.PARA.save_date(4:5));
                [year, ~,~] = datevec(tile.FORCING.PARA.start_time);
                
                if tile.FORCING.PARA.start_time <= datenum(year, save_month, save_day)
                    ground.TEMP.next_store_time = datenum(year + ground.PARA.start_store_year - 1, save_month, save_day) + ground.PARA.store_interval;
                    ground.TEMP.next_save_time = datenum(year + ground.PARA.start_store_year, save_month, save_day);
                    ground.TEMP.next_flip_time = datenum(year + ground.PARA.number_of_model_years, save_month, save_day);
                else
                    ground.TEMP.next_store_time = datenum(year + ground.PARA.start_store_year, save_month, save_day) + ground.PARA.store_interval;
                    ground.TEMP.next_save_time = datenum(year + ground.PARA.start_store_year + 1, save_month, save_day);
                    ground.TEMP.next_flip_time = datenum(year + ground.PARA.number_of_model_years +1, save_month, save_day);
                end
% 
%             if tile.FORCING.PARA.start_time <= datenum([ground.PARA.save_date datestr(tile.FORCING.PARA.start_time, 'yyyy')], 'dd.mm.yyyy')
%                 ground.TEMP.next_store_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time, 'yyyy')) + ground.PARA.start_store_year - 1), ], 'dd.mm.yyyy') + ground.PARA.store_interval;
%                 ground.TEMP.next_save_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time, 'yyyy')) + ground.PARA.start_store_year)], 'dd.mm.yyyy');
%                 ground.TEMP.next_flip_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time, 'yyyy')) + ground.PARA.number_of_model_years)], 'dd.mm.yyyy');
%             else
%                 ground.TEMP.next_store_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time, 'yyyy')) + ground.PARA.start_store_year)], 'dd.mm.yyyy') + ground.PARA.store_interval;
%                 ground.TEMP.next_save_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time, 'yyyy')) + ground.PARA.start_store_year + 1)], 'dd.mm.yyyy');
%                 ground.TEMP.next_flip_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time, 'yyyy')) + ground.PARA.number_of_model_years + 1)], 'dd.mm.yyyy');
%             end
            
            ground.TEMP.year_count = 1;
            ground.STORE.stratigraphy =[];
        end
        
        %----------------------
        
        function timestep = get_timestep_flip_flop(ground, tile)
                            
            timestep = (ground.TEMP.next_store_time - tile.t) .* ground.CONST.day_sec;

        end
        
        function timestep = get_timestep_flip_flop2(ground, tile)
                            
            timestep = (ground.TEMP.next_read_time - tile.t) .* ground.CONST.day_sec;

        end
        
        %----check_trigger
        %at next_store_time
        function ground = store_flip_flop(ground, tile)
            if tile.t == ground.TEMP.next_store_time
                copy_of_ground = GROUND_store_flip_flop_singleClass();
                copy_of_ground.STATVAR = ground.STATVAR;
                copy_of_ground.STATVAR.timestamp = tile.t;
                %copy_of_ground.PARA = ground.PARA;

                ground.STORE.stratigraphy{1,size(ground.STORE.stratigraphy,2)+1} = copy_of_ground;
                
                ground.TEMP.next_store_time = tile.t + ground.PARA.store_interval;
            end
        end        
        
        
        %at next_save_time
        function ground = write_store_flip_flop(ground, tile)
            if tile.t >= ground.TEMP.next_save_time
                run_name = tile.PARA.run_name;
                result_path = tile.PARA.result_path;
                if ~(exist([result_path run_name])==7)
                    mkdir([result_path run_name])
                end
                store = ground.STORE;
                save([result_path run_name '/' run_name '_STORE_' num2str(ground.TEMP.year_count) '.mat'], 'store')
                
                ground.STORE.stratigraphy = [];
                ground.TEMP.year_count = ground.TEMP.year_count + 1;
                
                save_day = str2num(ground.PARA.save_date(1:2));
                save_month = str2num(ground.PARA.save_date(4:5));
                [year, ~,~] = datevec(ground.TEMP.next_save_time);

                ground.TEMP.next_save_time = datenum(year + 1, save_month, save_day);
               % ground.TEMP.next_save_time = datenum([ground.PARA.save_date num2str(str2num(datestr(ground.TEMP.next_save_time, 'yyyy')) + 1)], 'dd.mm.yyyy');
            end
            
            
        end
        
        %at next_flip_time
        function ground = switch2read_STATVAR_flip_flop(ground, tile)
            %move "normal" class to STORE
            %change IA class downwards
             if tile.t >= ground.TEMP.next_flip_time
                 save_day = str2num(ground.PARA.save_date(1:2));
                 save_month = str2num(ground.PARA.save_date(4:5));
                 [year, ~,~] = datevec(tile.t);
                 ground.TEMP.next_flip_time = datenum(year + ground.PARA.number_of_read_years, save_month, save_day);
                 
                 %ground.TEMP.next_flip_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.t, 'yyyy')) + ground.PARA.number_of_read_years)], 'dd.mm.yyyy');
                 new_ground = GROUND_store_flip_flop_singleClass();
                 %new_ground.TEMP.next_flip_time = ground.TEMP.next_flip_time;
                 new_ground.TEMP = ground.TEMP;
                 new_ground.PARA = ground.PARA;
                 new_ground.STATVAR = ground.STATVAR;
                 new_ground.STATVAR.timestamp = tile.t;
                 new_ground.CONST = ground.CONST;
                 %new_ground.PARA.model_class = class(ground);
                 new_ground = finalize_init(new_ground, tile);
                 new_ground.NEXT = ground.NEXT;
                 new_ground.PREVIOUS = tile.TOP;  %becomes the first class in the STRATIGRAPHY
                 
                 new_ground.MODEL_CLASS = copy(ground);
                 new_ground.MODEL_CLASS.NEXT = [];
                 new_ground.MODEL_CLASS.PREVIOUS = [];
                 new_ground.MODEL_CLASS.IA_NEXT = [];
                 new_ground.MODEL_CLASS.IA_PREVIOUS = [];
                 
                 ground = new_ground;
                 ground.PREVIOUS.NEXT = ground;
                 ground.NEXT.PREVIOUS = ground;
                 ground.IA_PREVIOUS = [];
                 if ~isequal(ground.NEXT, tile.BOTTOM)
                     ground.IA_NEXT = get_IA_class(class(ground), class(ground.NEXT));
                     ground.NEXT.IA_PREVIOUS = ground.IA_NEXT;
                     ground.IA_NEXT.PREVIOUS = ground;
                     ground.IA_NEXT.NEXT = ground.NEXT;
                 end
                 tile.LATERAL.IA_TIME = ground.TEMP.next_flip_time + tile.LATERAL.IA_TIME_INCREMENT; %no lateral interactions during READ phase
            end
            
        end
        
        
        function ground = switch2read_STATVAR_flip_flop_BGC(ground, tile)
            %move "normal" class to STORE
            %change IA class downwards
             if tile.t >= ground.TEMP.next_flip_time
                 save_day = str2num(ground.PARA.save_date(1:2));
                 save_month = str2num(ground.PARA.save_date(4:5));
                 [year, ~,~] = datevec(tile.t);
                 ground.TEMP.next_flip_time = datenum(year + ground.PARA.number_of_read_years, save_month, save_day);
                 %ground.TEMP.next_flip_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.t, 'yyyy')) + ground.PARA.number_of_read_years)], 'dd.mm.yyyy');
                 
                 new_ground = GROUND_store_flip_flop_singleClass_BGC();
                 
                 new_ground.BGC = ground.BGC;
                 new_ground.BGC.PARENT = new_ground;
                 new_ground.IA_BGC = IA_BGC_read_statvar_from_out();
                 new_ground.IA_BGC.BGC = new_ground.BGC;
                 new_ground.IA_BGC.GROUND = new_ground;
                 new_ground.BGC.IA_BGC = new_ground.IA_BGC;
            
                 new_ground.TEMP = ground.TEMP;
                 new_ground.PARA = ground.PARA;
                 new_ground.STATVAR = ground.STATVAR;
                 new_ground.STATVAR.timestamp = tile.t;
                 new_ground.CONST = ground.CONST;
                 %new_ground.PARA.model_class = class(ground);
                 new_ground = finalize_init(new_ground, tile);
                 new_ground.NEXT = ground.NEXT;
                 new_ground.PREVIOUS = tile.TOP;  %becomes the first class in the STRATIGRAPHY
                 
                 new_ground.MODEL_CLASS = copy(ground);
                 new_ground.MODEL_CLASS.NEXT = [];
                 new_ground.MODEL_CLASS.PREVIOUS = [];
                 new_ground.MODEL_CLASS.IA_NEXT = [];
                 new_ground.MODEL_CLASS.IA_PREVIOUS = [];
                 new_ground.MODEL_CLASS.BGC = [];
                 new_ground.MODEL_CLASS.IA_BGC.BGC = []; %leave IA_BGC in place, but cut all dependnecies
                 new_ground.MODEL_CLASS.IA_BGC.GROUND = [];
                 
                 ground = new_ground;
                 ground.PREVIOUS.NEXT = ground;
                 ground.NEXT.PREVIOUS = ground;
                 ground.IA_PREVIOUS = [];
                 if ~isequal(ground.NEXT, tile.BOTTOM)
                     ground.IA_NEXT = get_IA_class(class(ground), class(ground.NEXT));
                     ground.NEXT.IA_PREVIOUS = ground.IA_NEXT;
                     ground.IA_NEXT.PREVIOUS = ground;
                     ground.IA_NEXT.NEXT = ground.NEXT;
                 end
                 
                 tile.LATERAL.IA_TIME = ground.TEMP.next_flip_time + tile.LATERAL.IA_TIME_INCREMENT; %no lateral interactions during READ phase
             end
            
            
        end
        
        %at next_flip_time
        function ground = switch2model_flip_flop(ground, tile)
            %take "normal" class from STORE and do nothing, STATVAR is still fine
            %change IA_CLASS downwards
            if tile.t >= ground.TEMP.next_flip_time
                
                new_ground = ground.MODEL_CLASS;
                
                new_ground.NEXT = ground.NEXT;
                new_ground.PREVIOUS = tile.TOP;  %becomes the first class in the STRATIGRAPHY
                ground = new_ground;
                ground.PREVIOUS.NEXT = ground;
                ground.NEXT.PREVIOUS = ground;
                ground.IA_PREVIOUS = [];
                if ~isequal(ground.NEXT, tile.BOTTOM)
                    ground.IA_NEXT = get_IA_class(class(ground), class(ground.NEXT));
                    ground.NEXT.IA_PREVIOUS = ground.IA_NEXT;
                    ground.IA_NEXT.PREVIOUS = ground;
                    ground.IA_NEXT.NEXT = ground.NEXT;
                end
                
                ground.TEMP.year_count = 1;
                ground.STORE.stratigraphy =[];
                
                save_day = str2num(ground.PARA.save_date(1:2));
                save_month = str2num(ground.PARA.save_date(4:5));
                [year, ~,~] = datevec(tile.t);
                ground.TEMP.next_store_time = max(datenum(year + ground.PARA.start_store_year - 1, save_month, save_day) + ground.PARA.store_interval, tile.t + ground.PARA.store_interval);
                ground.TEMP.next_save_time = datenum(year + ground.PARA.start_store_year, save_month, save_day);
                ground.TEMP.next_flip_time = datenum(year + ground.PARA.number_of_model_years, save_month, save_day);
                
%                 ground.TEMP.next_store_time = max(datenum([ground.PARA.save_date num2str(str2num(datestr(tile.t, 'yyyy')) + ground.PARA.start_store_year - 1), ], 'dd.mm.yyyy') + ground.PARA.store_interval, tile.t + ground.PARA.store_interval);
%                 ground.TEMP.next_save_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.t, 'yyyy')) + ground.PARA.start_store_year)], 'dd.mm.yyyy');
%                 ground.TEMP.next_flip_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.t, 'yyyy')) + ground.PARA.number_of_model_years)], 'dd.mm.yyyy');
            end
            
            
        end
        
        %at next_flip_time
        function ground = switch2model_flip_flop_BGC(ground, tile)
            %take "normal" class from STORE and do nothing, STATVAR is still fine
            %change IA_CLASS downwards
            if tile.t >= ground.TEMP.next_flip_time
                
                new_ground = ground.MODEL_CLASS;
                
                new_ground.BGC = ground.BGC; %switch BGC anmd reconnect IA_BGC class
                new_ground.BGC.PARENT = new_ground;
                new_ground.IA_BGC.GROUND = new_ground;
                new_ground.IA_BGC.BGC = new_ground.BGC;
                new_ground.BGC.IA_BGC = new_ground.IA_BGC;
                %new_ground.BGC.TEMP.regrid_flag = 1;
               
                 new_ground.NEXT = ground.NEXT;
                 new_ground.PREVIOUS = tile.TOP;  %becomes the first class in the STRATIGRAPHY
                 ground = new_ground;
                 ground.PREVIOUS.NEXT = ground;
                 ground.NEXT.PREVIOUS = ground;
                 ground.IA_PREVIOUS = [];
                 if ~isequal(ground.NEXT, tile.BOTTOM)
                     ground.IA_NEXT = get_IA_class(class(ground), class(ground.NEXT));
                     ground.NEXT.IA_PREVIOUS = ground.IA_NEXT;
                     ground.IA_NEXT.PREVIOUS = ground;
                     ground.IA_NEXT.NEXT = ground.NEXT;
                 end
                 
                 add_new_BGC_cells_from_READ(ground.IA_BGC, tile);
                 
                 ground.TEMP.year_count = 1;
                 ground.STORE.stratigraphy =[];
                                  
                 save_day = str2num(ground.PARA.save_date(1:2));
                 save_month = str2num(ground.PARA.save_date(4:5));
                 [year, ~,~] = datevec(tile.t);
                 ground.TEMP.next_store_time = max(datenum(year + ground.PARA.start_store_year - 1, save_month, save_day) + ground.PARA.store_interval, tile.t + ground.PARA.store_interval);
                 ground.TEMP.next_save_time = datenum(year + ground.PARA.start_store_year, save_month, save_day);
                 ground.TEMP.next_flip_time = datenum(year + ground.PARA.number_of_model_years, save_month, save_day);
%                  
%                  ground.TEMP.next_store_time = max(datenum([ground.PARA.save_date num2str(str2num(datestr(tile.t, 'yyyy')) + ground.PARA.start_store_year - 1), ], 'dd.mm.yyyy') + ground.PARA.store_interval, tile.t + ground.PARA.store_interval);
%                  ground.TEMP.next_save_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.t, 'yyyy')) + ground.PARA.start_store_year)], 'dd.mm.yyyy');
%                  ground.TEMP.next_flip_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.t, 'yyyy')) + ground.PARA.number_of_model_years)], 'dd.mm.yyyy');
            end
        end
        
    end
end

