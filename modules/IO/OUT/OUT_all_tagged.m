%========================================================================
% CryoGrid OUT class defining storage format of the output 
% OUT_all_tagged stores identical copies of all GROUND classses (including 
% STATVAR, TEMP, PARA) in the stratigraphy for each output timestep, while 
% lateral classes are not stored.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
%
% This class takes an additional PARA parameter 'tag' which is added to the
% output filename, such that the final output filenames are:
% [run_name '_' tag '_' date '.mat'].
% If 'tag' is not specified in PARA (or blank), it will fall back to the
% standard output naming.
% If save interval is not specified, it will only store results and the end
% of the entire model run.
%
% T. Ingeman-Nielsen, S. Westermann, J. Scheer, December 2020
%========================================================================


classdef OUT_all_tagged < matlab.mixin.Copyable
 

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
        
        function out = initialize_excel(out)
            
        end
        
        
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
        
        %---------------time integration-------------
		
% 		function out = store_OUT(out, t, TOP, BOTTOM, forcing, run_number, timestep, result_path)
            
        function out = store_OUT(out, tile)           
            
             t = tile.t;
             TOP = tile.TOP; 
             BOTTOM = tile.BOTTOM;
             forcing = tile.FORCING;
             run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
             result_path = tile.PARA.result_path;            
             timestep = tile.timestep;
             out_tag = out.PARA.tag;

            
            
            if t==out.OUTPUT_TIME
                %if id == 1
                disp([datestr(t)])
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
                    res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_NEXT=[];  %cut all dependencies
                    if isprop(res, 'CHILD')
                        res.CHILD = [];
                        res.IA_CHILD =[];
                    end
                    result=[result; {res}];
                    CURRENT = CURRENT.NEXT;
                end
                out.STRATIGRAPHY{1,size(out.STRATIGRAPHY,2)+1} = result;
                
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
                    out.TIMESTAMP=[];
                    out.MISC=[];
                    if ~isnan(out.PARA.save_interval)
                        out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                    end
                end
            end
        end

        function xls_out = write_excel(out)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			
            xls_out = {'OUT','index',NaN,NaN;'OUT_all',1,NaN,NaN;'output_timestep',0.250000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'OUT_END',NaN,NaN,NaN};
        end
        
    end
end