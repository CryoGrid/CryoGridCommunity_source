%========================================================================
% CryoGrid OUT class defining storage format of the output 
% OUT_all_lateral stores identical copies of all GROUND classses (including STATVAR, TEMP, PARA) in the
% stratigraphy for each output timestep and copies of all LATERAL classes.
% Other than that, it is identical to OUT_all.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================


classdef OUT_all_lateral_vegetation < matlab.mixin.Copyable

    properties
		out_index
        STRATIGRAPHY
        LATERAL
        VEGETATION
        FORCING
        TIMESTAMP
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
		CONST
	end
    
    
    methods
		

% 
%         function out = initialize_excel(out)
%             
%         end
        
        function out = provide_PARA(out)         
            % INITIALIZE_PARA  Initializes PARA structure.

            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
        end
        
        function out = provide_CONST(out)

        end
        
        function out = provide_STATVAR(out)

        end
		
% 	    function out = initialize_TEMP(out)
%             % INITIALIZE_TEMP  Initializes TEMP structure.
%             
%             out.TEMP = struct();
%         end	
	
% 		function out = populate_PARA(out, pprovider)
%             % POPULATE_PARA  Updates the PARA structure with values from pprovider.
%             %
%             %   ARGUMENTS:
%             %   pprovider:  instance of PARAMETER_PROVIDER class
%             
%             out.PARA = pprovider.populate_struct(out.PARA, 'OUT', mfilename('class'), out.out_index);
%         end
		
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
        end
        
        %-------time integration----------------
		
		%function out = store_OUT(out, t, TOP, BOTTOM, forcing, run_number, timestep, result_path)
            
        function out = store_OUT(out, tile)
            
             t = tile.t;
             TOP = tile.TOP; 
             BOTTOM = tile.BOTTOM;
             forcing = tile.FORCING;
             %run_number = tile.RUN_NUMBER;
             run_name = tile.PARA.run_name;
             result_path = tile.PARA.result_path;
             timestep = tile.timestep;
             out_tag = out.PARA.tag;
             
            
            if t==out.OUTPUT_TIME
                %if id == 1
                disp([datestr(t)])
                %end
                %labBarrier
                out.TIMESTAMP=[out.TIMESTAMP t];
                
                out.FORCING{1,size(out.FORCING, 2)+1} = tile.FORCING.TEMP;
%                 OUT.FORCING = OUT.FORCING ;
                
                
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
                        res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_PREVIOUS=[];  res.PARENT = []; %cut all dependencies
                        result=[result; {res}];
                    end
                    res = copy(CURRENT);
                    if isprop(res, 'VEGETATION')
                        res.VEGETATION =[];  %remove look-up tables, runs out of memory otherwise
                    end
                    if isprop(res, 'LUT')
                        res.LUT =[];  %remove look-up tables, runs out of memory otherwise
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
                
                %lateral, read only STATVAR and PARA---
                result={};
                ia_classes=TOP.LATERAL.IA_CLASSES;
                for i=1:size(ia_classes,1)
                    res = copy(ia_classes{i,1});
                    vars = fieldnames(res);
                    for j=1:size(vars,1)
                        if ~strcmp(vars{j,1}, 'PARA') && ~strcmp(vars{j,1}, 'STATVAR')
                           res.(vars{j,1}) = []; 
                        end
                    end
                    result=[result; {res}];                    
                end
                out.LATERAL{1,size(out.LATERAL, 2)+1} = result;
                
                if isprop(TOP.NEXT, 'VEGETATION')
                    out.VEGETATION{1,size(out.VEGETATION, 2)+1} = copy(TOP.NEXT.VEGETATION);
                end
                    
                
                %---
                
                
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
                if t==out.SAVE_TIME 
                   if ~(exist([result_path run_name])==7)
                       mkdir([result_path run_name])
                   end
                   if isempty(out_tag) || all(isnan(out_tag))
                       save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '.mat'], 'out')
                   else
                       save([result_path run_name '/' run_name '_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'out')
                   end
                   
                   out.STRATIGRAPHY=[];
                   out.VEGETATION =[];
                   out.FORCING = [];
                   out.LATERAL=[];
                   out.TIMESTAMP=[];
                   out.MISC=[];
                   out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end
            end
        end
        
        %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'OUT';
           
            out.PARA.default_value.output_timestep = {0.25};
            out.PARA.comment.output_timestep = {'timestep of output [days]'};

            out.PARA.default_value.save_date = {'01.09.'};
            out.PARA.comment.save_date = {'date (dd.mm.) when output file is written'};
            
            out.PARA.default_value.save_interval = {1};
            out.PARA.comment.save_interval = {'interval of output files [years]'};
            
            out.PARA.default_value.tag = {''};
            out.PARA.comment.tag = {'additional tag added to file name'};
        end

%         function xls_out = write_excel(out)
% 			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
% 			
%             xls_out = {'OUT','index',NaN,NaN;'OUT_all',1,NaN,NaN;'output_timestep',0.250000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'OUT_END',NaN,NaN,NaN};
%         end
         
        
% 		% ==========================================
%         % DEPRECATED METHODS
%         % to be deleted when new implementation
%         % is validated for backwards compatibility
%         % ==========================================
% 		
% 
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