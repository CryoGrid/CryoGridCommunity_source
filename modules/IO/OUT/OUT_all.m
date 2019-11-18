% forcing data


classdef OUT_all
    properties
        STRATIGRAPHY
        TIMESTAMP
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
    end
    
    
    methods
        
        function res = initialize_OUT_all(run_info) %discontinued
            res.PARA.output_timestep = 1/4;
            res.PARA.save_date = '01.04.';
            res.PARA.save_interval = 0.1;
            res.OUTPUT_TIME = run_info.START_TIME + res.PARA.output_timestep;
            if isempty (res.PARA.save_interval)
                res.SAVE_TIME = run_info.END_TIME;
            else
                res.SAVE_TIME = min(run_info.END_TIME,  datenum([res.PARA.save_date num2str(str2num(datestr(run_info.START_TIME,'yyyy')) + res.PARA.save_interval)], 'dd.mm.yyyy'));
            end
        end
        
        function out = provide_variables(out)
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
        end
        
        function out = initalize_from_file(out, section)
           variables = fieldnames(out.PARA);
           for i=1:size(variables,1)
               for j=1:size(section,1)
                  if strcmp(variables{i,1}, section{j,1})
                      out.PARA.(variables{i,1}) = section{j,2};
                  end
               end
           end
        end
        
        function out = complete_init_out(out, forcing)
            out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));

%                 out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(datenum(forcing.PARA.start_time,'dd.mm.yyyy'),'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
        end
            
        
        function out = store_OUT(out, t, TOP_CLASS, BOTTOM, forcing, run_number)
            if t==out.OUTPUT_TIME
                  if TOP_CLASS == GROUND_vegetation_snow
%                     disp(datestr(t))
                    
                    out.TIMESTAMP=[out.TIMESTAMP t];
                    
                    %                 out.STRATIGRAPHY{1,size(out.STRATIGRAPHY,2)+1} = copy(TOP_CLASS);  %append new stratigraphy, should be made more sophisticated by not adding instaneous values, but averaging/accumulating variables
                    % tveg (sunny) + tground
                    out.STRATIGRAPHY{size(out.STRATIGRAPHY,1)+1,2} = [TOP_CLASS.STATVAR.vegetation.mlcanopyinst.tveg(:,:,1)'; TOP_CLASS.NEXT.STATVAR.T(:,:)];
                    % lhflx, shflx, storageflx, net radiation, photosynthesis, fraction sunny/shaded
                    out.STRATIGRAPHY{size(out.STRATIGRAPHY,1),3} = TOP_CLASS.STATVAR.vegetation.mlcanopyinst.st_prof';
                    out.STRATIGRAPHY{size(out.STRATIGRAPHY,1),4} = TOP_CLASS.STATVAR.vegetation.mlcanopyinst.sh_prof';
                    out.STRATIGRAPHY{size(out.STRATIGRAPHY,1),5} = TOP_CLASS.STATVAR.vegetation.mlcanopyinst.lh_prof';
                    out.STRATIGRAPHY{size(out.STRATIGRAPHY,1),6} = TOP_CLASS.STATVAR.vegetation.mlcanopyinst.rn_prof';
                    out.STRATIGRAPHY{size(out.STRATIGRAPHY,1),7} = TOP_CLASS.STATVAR.vegetation.mlcanopyinst.an(:,:,1)';
                    out.STRATIGRAPHY{size(out.STRATIGRAPHY,1),8} = TOP_CLASS.STATVAR.vegetation.flux.fracsun';
                    out.STRATIGRAPHY{size(out.STRATIGRAPHY,1),9} = TOP_CLASS.STATVAR.vegetation.flux.fracsha';
                    % wind, precipitation
                    %                 out.STRATIGRAPHY{size(out.STRATIGRAPHY,1),10} = TOP_CLASS.STATVAR.vegetation.mlcanopyinst.wind';
                    
                    % tveg
                    %                 out.STRATIGRAPHY{size(out.STRATIGRAPHY,1)+1,2} = [TOP_CLASS.STATVAR.vegetation.mlcanopyinst.tveg(:,:,1)'; TOP_CLASS.NEXT.STATVAR.T(:,:)];
                    
                    CURRENT = TOP_CLASS;
                    if isprop(CURRENT, 'IA_CHILD') && ~isempty(CURRENT.IA_CHILD)
                        out.MISC=[out.MISC [CURRENT.IA_CHILD.IA_CHILD_SNOW.STATVAR.T(1,1); CURRENT.IA_CHILD.IA_CHILD_SNOW.STATVAR.layerThick(1,1)]];
                    else
                        out.MISC=[out.MISC [NaN; NaN]];
                    end
                    result={};
                    while ~isequal(CURRENT, BOTTOM)
                        if isprop(CURRENT, 'IA_CHILD') && ~isempty(CURRENT.IA_CHILD) && CURRENT.IA_CHILD.STATUS ==1
                            res=copy(CURRENT.IA_CHILD.IA_CHILD_SNOW);
                            res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_NEXT=[];  %cut all dependencies
                            result=[result; {res}];
                        end
                        res = copy(CURRENT);
                        res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_NEXT=[];  %cut all dependencies
                        if isprop(res, 'IA_CHILD')
                            res.IA_CHILD =[];
                        end
                        result=[result; {res}];
                        CURRENT = CURRENT.NEXT;
                    end
                    out.STRATIGRAPHY{size(out.STRATIGRAPHY,1),1} = result;
                    
                    out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
                    if t==out.SAVE_TIME
                        
                        save(['results/' run_number '/' run_number '_' datestr(t,'yyyy') '.mat'], 'out')
                        out.STRATIGRAPHY=[];
                        out.TIMESTAMP=[];
                        out.MISC=[];
                        out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                    end
                else
%                     disp(datestr(t))
                    
                    out.TIMESTAMP=[out.TIMESTAMP t];
                    
                    %out.STRATIGRAPHY{1,size(out.STRATIGRAPHY,2)+1} = copy(TOP_CLASS);  %append new stratigraphy, should be made more sophisticated by not adding instaneous values, but averaging/accumulating variables
                    %out.STRATIGRAPHY{1,size(out.STRATIGRAPHY,2)+1} = [TOP_CLASS.STATVAR.T; TOP_CLASS.NEXT.STATVAR.T];
                    CURRENT =TOP_CLASS;
                    if isprop(CURRENT, 'IA_CHILD') && ~isempty(CURRENT.IA_CHILD)
                        out.MISC=[out.MISC [CURRENT.IA_CHILD.IA_CHILD_SNOW.STATVAR.T(1,1); CURRENT.IA_CHILD.IA_CHILD_SNOW.STATVAR.layerThick(1,1)]];
                    else
                        out.MISC=[out.MISC [NaN; NaN]];
                    end
                    result={};
                    while ~isequal(CURRENT, BOTTOM)
                        if isprop(CURRENT, 'IA_CHILD') && ~isempty(CURRENT.IA_CHILD) && CURRENT.IA_CHILD.STATUS ==1
                            res=copy(CURRENT.IA_CHILD.IA_CHILD_SNOW);
                            res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_NEXT=[];  %cut all dependencies
                            result=[result; {res}];
                        end
                        res = copy(CURRENT);
                        res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_NEXT=[];  %cut all dependencies
                        if isprop(res, 'IA_CHILD')
                            res.IA_CHILD =[];
                        end
                        result=[result; {res}];
                        CURRENT = CURRENT.NEXT;
                    end
                    out.STRATIGRAPHY{1,size(out.STRATIGRAPHY,2)+1} = result;
                    
                    out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
                    if t==out.SAVE_TIME
                        
                        save(['results/' run_number '/' run_number '_' datestr(t,'yyyy') '.mat'], 'out')
                        out.STRATIGRAPHY=[];
                        out.TIMESTAMP=[];
                        out.MISC=[];
                        out.SAVE_TIME = forcing.PARA.end_time; %min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                    end
                    
                end
            end
        end
        
    end
end