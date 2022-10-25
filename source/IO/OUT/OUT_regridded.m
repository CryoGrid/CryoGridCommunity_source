%========================================================================
% CryoGrid OUT class OUT_all
% CryoGrid OUT class defining storage format of the output 
% OUT_all stores identical copies of all GROUND classses (including STATVAR, TEMP, PARA) in the
% stratigraphy for each output timestep, while lateral classes are not stored.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, June 2021
%========================================================================


classdef OUT_regridded < matlab.mixin.Copyable
 

    properties
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
        CONST
        result
        timestamp
        
    end
    
    
    methods
        
        %initialization
        
        function out = provide_PARA(out)         

            out.PARA.variables = [];
            out.PARA.upper_elevation = [];
            out.PARA.lower_elevation = [];
            out.PARA.target_grid_size = [];
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

            forcing = tile.FORCING;
            
            % Set the next (first) output time. This is the next (first) time output
            % is collected (in memory) for later storage to disk.
            out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            %make struct "result" and initialize all variables defined by
            %the user as empty arrays
            out.timestamp = [];
            for i =  1:size(out.PARA.variables,1)
               out.result.(out.PARA.variables{i,1}) = []; 
            end
            out.result.depths = [];
            out.result.class_number = [];
            
            % Set the next (first) save time. This is the next (first) time all the
            % collected output is saved to disk.
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
                        
        end
        
        %---------------time integration-------------
        
%         function out = store_OUT(out, t, TOP, BOTTOM, forcing, run_number, timestep, result_path)
            
        function out = store_OUT(out, tile)           
            
            t = tile.t;
            TOP = tile.TOP; 
            BOTTOM = tile.BOTTOM;
            forcing = tile.FORCING;
            run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
            result_path = tile.PARA.result_path;            
            timestep = tile.timestep;
            out_tag = out.PARA.tag;
            
            if t>=out.OUTPUT_TIME
                % It is time to collect output
                % Store the current state of the model in the out structure.

                disp([datestr(t)])                
                out.timestamp = [out.timestamp t];
                
                CURRENT =TOP.NEXT;
                if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
                    out.MISC= [CURRENT.CHILD.STATVAR.T(1,1); CURRENT.CHILD.STATVAR.layerThick(1,1)]; 
                else
                    out.MISC = [NaN; NaN];
                end
                result={};
                while ~isequal(CURRENT, BOTTOM)
                    if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
                        res=copy(CURRENT.CHILD);
                        res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_PREVIOUS=[];  res.PARENT = []; %cut all dependencies
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
                %result contains the information needed for regridding 
                %regrid results and append to out.result
                out = regrid_out(out, result);
                
                % Set the next OUTPUT_TIME
                out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                
                if t>=out.SAVE_TIME
                    % It is time to save all the collected model output to disk
                     
                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end
                    CG_out = out.result;
                    CG_out.timestamp = out.timestamp;
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    else
                        save([result_path run_name '/' run_name '_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    end
                    
                    % Clear the out structure
                    out.MISC=[];
                    %make struct "result" and initialize all variables defined by
                    %the user as empty arrays
                    out.timestamp = [];
                    for i =  1:size(out.PARA.variables,1)
                        out.result.(out.PARA.variables{i,1}) = [];
                    end
                    out.result.depths = [];
                    out.result.class_number = [];
                    
                    if ~isnan(out.PARA.save_interval)
                        % If save_interval is defined, uptate SAVE_TIME for next save opertion 
                        out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        % If save_interval is not defined, we will save at the very end of the model run
                        % and thus do not need to update SAVE_TIME (update would fail because save_interval is nan)
					end
                end
            end
        end
        
        %taken from read_display_out
        function out = regrid_out(out, result)
            new_grid = [out.PARA.upper_elevation:-out.PARA.target_grid_size:out.PARA.lower_elevation]';
            threshold = out.PARA.target_grid_size/10;
            
            variableList = fieldnames(out.result);
            numberOfVariables = size(variableList,1);
            
            altitudeLowestCell = result{end,1}.STATVAR.lowerPos;

            %read out and accumulate over all classes
            
            layerThick=[];
            area=[];
            
            for j=1:size(result,1)
                layerThick=[layerThick; result{j,1}.STATVAR.layerThick];
                area=[area; result{j,1}.STATVAR.area];
            end
            layerThick_temp = layerThick;
            layerThick = zeros(size(layerThick,1).*2, 1).* NaN;
            layerThick(1:2:size(layerThick,1),1) = threshold;
            layerThick(2:2:size(layerThick,1),1) = layerThick_temp - threshold;
            
            area_temp = area;
            area=[area; area].* NaN;
            area(1:2:size(layerThick,1),1) = area_temp;
            area(2:2:size(layerThick,1),1) = area_temp;
               
            temp=repmat(NaN, size(layerThick_temp,1), numberOfVariables);
            pos=1;
            for j = 1:size(result,1)
                fieldLength = size(result{j,1}.STATVAR.layerThick,1);
                for k=1:numberOfVariables-1
                    if any(strcmp(fieldnames(result{j,1}.STATVAR), variableList{k,1}))
                        temp(pos:pos+fieldLength-1,k) = result{j,1}.STATVAR.(variableList{k,1});
                    end
                end
                temp(pos:pos+fieldLength-1,numberOfVariables) = zeros(fieldLength,1) + size(result,1)+1-j; %assigna class number starting with 1 from the bottom
                pos = pos+fieldLength;
            end
            
            %compute targate variables
            for k=1:numberOfVariables
                if strcmp(variableList{k,1}, 'saltConc')
                    pos_waterIce = find(strcmp(variableList, 'waterIce'));
                    temp(:,k) = temp(:,k)./ (temp(:,pos_waterIce) ./layerThick_temp./area_temp);  %divide by total water content
                end
            end
            
            for k=1:numberOfVariables
                if strcmp(variableList{k,1}, 'water') || strcmp(variableList{k,1}, 'ice') || strcmp(variableList{k,1}, 'waterIce') || strcmp(variableList{k,1}, 'XwaterIce') || strcmp(variableList{k,1}, 'Xwater') || strcmp(variableList{k,1}, 'Xice') || strcmp(variableList{k,1}, 'saltConc')
                    temp(:,k) = temp(:,k)./layerThick_temp./area_temp;
                end
            end
            
            temp_temp = temp;
            temp=[temp; temp].* NaN;
            temp(1:2:size(layerThick,1),:) = temp_temp;
            temp(2:2:size(layerThick,1),:) = temp_temp;
            
            %interpolate to new grid
            depths = cumsum(layerThick);
            depths = depths - threshold/2;
            depths(1) = 0;
            depths = -(depths-depths(end,1));
            depths = depths + altitudeLowestCell;
            
            for k=1:numberOfVariables
                out.result.(variableList{k,1}) = [out.result.(variableList{k,1}) interp1(depths, temp(:,k), new_grid, 'nearest')];
            end
            out.result.depths = new_grid;
            
        end
        

        %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'OUT';
            
            out.PARA.comment.variables = {'select output variables: T, water, ice, waterIce, XwaterIce, Xwater, Xice, saltConc are supported'};
            out.PARA.options.variables.name = 'H_LIST';
            out.PARA.options.variables.entries_x = {'T' 'water' 'ice' 'waterIce'};
      
            out.PARA.comment.upper_elevation = {'upper elevation of output domain, in m a.s.l.; must match altitude in TILE!'};
            
            out.PARA.comment.lower_elevation = {'lower elevation of output domain, in m a.s.l.; must match altitude in TILE!'};
            
            out.PARA.default_value.target_grid_size = {0.02};
            out.PARA.comment.target_grid_size = {'cell size of regular grid that values are interpolated to'};
           
            out.PARA.default_value.output_timestep = {0.25};
            out.PARA.comment.output_timestep = {'timestep of output [days]'};

            out.PARA.default_value.save_date = {'01.09.'};
            out.PARA.comment.save_date = {'date (dd.mm.) when output file is written'};
            
            out.PARA.default_value.save_interval = {1};
            out.PARA.comment.save_interval = {'interval of output files [years]'};
            
            out.PARA.default_value.tag = {''};
            out.PARA.comment.tag = {'additional tag added to file name'};
        end
        

    end
end