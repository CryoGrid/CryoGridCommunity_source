%========================================================================
% CryoGrid OUT class defining storage format of the output 
% OUT_all stores identical copies of all GROUND classses (including STATVAR, TEMP, PARA) in the
% stratigraphy for each output timestep, while lateral classes are not stored.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================


classdef OUT_all < matlab.mixin.Copyable
 

    properties
		out_index
        STRATIGRAPHY
        LATERAL
        TIMESTAMP
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
        CONST
		
	end
    
    
    methods
		
        %initialization

        
        
        function out = provide_PARA(out)         

            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
        end
        
        function out = provide_CONST(out)

        end
        
        function out = provide_STATVAR(out)

        end
		

		
		function out = finalize_init(out, tile)
			% FINALIZE_SETUP  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
			
			%	ARGUMENTS:
			%	forcing:	instance of FORCING class
			forcing = tile.FORCING;
            
			out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
            
            out.TEMP = struct();
            out.TEMP.count = 0;
            out.TEMP.time = 0;
            
        end
        
        %---------------time integration-------------
		
% 		function out = store_OUT(out, t, TOP, BOTTOM, forcing, run_number, timestep, result_path)
            
        function out = store_OUT(out, tile)           
            
             t = tile.t;
             TOP = tile.TOP; 
             BOTTOM = tile.BOTTOM;
             forcing = tile.FORCING;
             run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
             result_path = tile.PARA.result_path;            

             out.TEMP.count = out.TEMP.count + 1;
             out.TEMP.time = out.TEMP.time + tile.timestep;
            
%           --- Write output ---
            if t>=out.OUTPUT_TIME
                %if id == 1
                avg_timestep = out.TEMP.time/out.TEMP.count;
                disp([datestr(t,'dd-mmm-yyyy HH:MM') ' Average timestep: ' num2str(avg_timestep)])
                out.TEMP.count = 0;
                out.TEMP.time = 0;
                %end
                %labBarrier
                out.TIMESTAMP=[out.TIMESTAMP t];
                
                CURRENT =TOP.NEXT;
                if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
                    out.MISC=[out.MISC [CURRENT.CHILD.STATVAR.T(1,1); CURRENT.CHILD.STATVAR.layerThick(1,1)]]; 
                else
                    out.MISC=[out.MISC [NaN; NaN]];
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
                    if isprop(res, 'READ_OUT')
                        res.READ_OUT =[];  %remove look-up tables, runs out of memory otherwise
                    end
                    res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_PREVIOUS=[];  %cut all dependencies
                    if isprop(res, 'CHILD')
                        res.CHILD = [];
                        res.IA_CHILD =[];
                    end
                    result=[result; {res}];
                    CURRENT = CURRENT.NEXT;
                end
                out.STRATIGRAPHY{1,size(out.STRATIGRAPHY,2)+1} = result;
                
                out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);

%               --- Store output ---
                if t>=out.SAVE_TIME                     
                   if ~(exist([result_path run_name])==7)
                       mkdir([result_path run_name])
                   end
                   save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '.mat'], 'out')
                   out.STRATIGRAPHY=[];
                   out.TIMESTAMP=[];
                   out.MISC=[];
                   out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end
            end
        end

        function xls_out = write_excel(out)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			
            xls_out = {'OUT','index',NaN,NaN;'OUT_all',1,NaN,NaN;'output_timestep',0.250000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'OUT_END',NaN,NaN,NaN};
        end
         
        
%         		function out = initialize(out)
%             % INITIALIZE  Initializes all properties needed by the class.
% 			
% 			out.STRATIGRAPHY = [];
% 			out.TIMESTAMP = [];
% 			out.MISC = [];
% 			out.OUTPUT_TIME = [];
% 			out.SAVE_TIME = [];
%             out = out.initialize_PARA();
%             out = out.initialize_TEMP();
%         end

% 	    function out = initialize_TEMP(out)
%             % INITIALIZE_TEMP  Initializes TEMP structure.
%             
%             out.TEMP = struct();
%         end	
        
        
		
%         function out = populate_PARA(out, pprovider)
%             % POPULATE_PARA  Updates the PARA structure with values from pprovider.
%             %
%             %   ARGUMENTS:
%             %   pprovider:  instance of PARAMETER_PROVIDER class
%             
%             out.PARA = pprovider.populate_struct(out.PARA, 'OUT', mfilename('class'), out.out_index);
%         end
        
% 		% ==========================================
%         % DEPRECATED METHODS
%         % to be deleted when new implementation
%         % is validated for backwards compatibility
%         % ==========================================
% 		
% 		% function res = initialize_out_all(run_info) %discontinued
% %             res.PARA.output_timestep = 1/4;
% %             res.PARA.save_date = '01.09.';
% %             res.PARA.save_interval = 1;
% %             res.OUTPUT_TIME = run_info.START_TIME + res.PARA.output_timestep;
% %             if isempty (res.PARA.save_interval)
% %                 res.SAVE_TIME = run_info.END_TIME;
% %             else
% %                 res.SAVE_TIME = min(run_info.END_TIME,  datenum([res.PARA.save_date num2str(str2num(datestr(run_info.START_TIME,'yyyy')) + res.PARA.save_interval)], 'dd.mm.yyyy'));
% %             end
% %         end
% 		
%         function out = provide_variables(out)
% 			st = dbstack;
%             warning(['DEPRECATION WARNING: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Use PARAMETER_PROVIDER class to obtain parameter values.']);
%             out.PARA.output_timestep = [];
%             out.PARA.save_date = [];
%             out.PARA.save_interval = [];
%         end
%         
%         function out = initalize_from_file(out, section)
% 			st = dbstack;
% 			warning(['DEPRECATION WARNING: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Use PARAMETER_PROVIDER class to obtain parameter values.']);
% 					 
% 			variables = fieldnames(out.PARA);
% 			for i=1:size(variables,1)
% 				for j=1:size(section,1)
% 					if strcmp(variables{i,1}, section{j,1})
% 						out.PARA.(variables{i,1}) = section{j,2};
% 					end
% 				end
% 			end
% 		end
% 
% 		function out = initialize_from_ini(out, ini)
% 			% INITIALIZE_FROM_INI  Initializes the variables from output structure of the ini parser, and compares the
% 			%	names of the variables from the class to the ini structure.
% 			% 	If the variables from the class mismatch the variables from
% 			% 	the ini file, an error message is displayed.
% 			
% 			%	ARGUMENTS:
% 			%	ini:	output structure from the ini parser
% 					
% 			st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%                      'Code should be moved to new PARAMETER_PROVIDER class ',...
%                      'to streamline file access and the population of parameters.']);
% 					 
% 			ini_variables = fields(ini.OUT_all);
% 			out_variables = fieldnames(out.PARA);
% 			ismatch_class_ini_variables(ini_variables, out_variables) 
% 			for i=1:length(out_variables)
% 				for j=1:length(ini_variables)
% 					if strcmp(out_variables{i,1},ini_variables{i,1})
% 						out.PARA.(out_variables{i,1}) = ini.out_all.(ini_variables{i,1})
% 					end
% 				end
% 			end
% 		end
% 			
% 		% 
%         function out = complete_init_out(out, forcing)
% 			% COMPLETE_INIT_OUT  Completes the initialization of the member variables of the class by creating them based on out.PARA and forcing.PARA.
% 			
% 			%	ARGUMENTS:
% 			%	forcing:	instance of FORCING class
% 			
% 			st = dbstack;
%             warning(['DEPRECATION: Method ' st.name '() is deprecated and will be removed.' newline,...
%             'Parameter initialization should be finalized in the ' mfilename('class') '.finalize_setup() ']);
%             
%             out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
%             if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
%                 out.SAVE_TIME = forcing.PARA.end_time;
%             else
%                 out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
%             end
%         end 
        
    end
end