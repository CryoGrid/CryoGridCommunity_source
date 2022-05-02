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


classdef OUT_all_lateral_CG3style
    
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
        
    end
    
    
    methods
        
        %constructor
        function out = OUT_all_lateral_CG3style(varargin)               % Temporary definition, to allow old code to run
            %function out = OUT_all(index, pprovider, forcing)      % Definition to be used when old code is no longer supported
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
            
            out.out_index = index;
            out = out.initialize();
            out = out.populate_PARA(pprovider);
            out = out.finalize_setup(forcing);
        end
        
        %-------initialization--------------
        
        function out = initialize(out)
            % INITIALIZE  Initializes all properties needed by the class.
            
            out.STRATIGRAPHY = [];
            out.LATERAL=[];
            out.TIMESTAMP = [];
            out.MISC = [];
            out.OUTPUT_TIME = [];
            out.SAVE_TIME = [];
            out = out.initialize_PARA();
            out = out.initialize_TEMP();
        end
        
        function out = initialize_PARA(out)
            % INITIALIZE_PARA  Initializes PARA structure.
            
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
        end
        
        function out = initialize_TEMP(out)
            % INITIALIZE_TEMP  Initializes TEMP structure.
            
            out.TEMP = struct();
        end
        
        function out = populate_PARA(out, pprovider)
            % POPULATE_PARA  Updates the PARA structure with values from pprovider.
            %
            %   ARGUMENTS:
            %   pprovider:  instance of PARAMETER_PROVIDER class
            
            out.PARA = pprovider.populate_struct(out.PARA, 'OUT', mfilename('class'), out.out_index);
        end
        
        function out = finalize_setup(out, forcing)
            % FINALIZE_SETUP  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
            
            %	ARGUMENTS:
            %	forcing:	instance of FORCING class
            
            out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval)
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
        end
        
        %-------time integration----------------
        
        %function out = store_OUT(out, t, TOP, BOTTOM, forcing, run_number, timestep, result_path)
        
        function out = store_OUT(out, tile)
            
            t = tile.t;
            TOP = tile.TOP;
            BOTTOM = tile.BOTTOM;
            forcing = tile.FORCING;
            run_number = tile.RUN_NUMBER;
            timestep = tile.timestep;
            result_path = tile.RESULT_PATH;
            
            
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
                        res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_PREVIOUS=[];  res.PARENT = []; %cut all dependencies
                        result=[result; {res}];
                    end
                    res = copy(CURRENT);
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
                %---
                
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
                if t==out.SAVE_TIME
                    if ~(exist([result_path run_number])==7)
                        mkdir([result_path run_number])
                    end
                    out2 = usableOUT(out);
                    save([result_path run_number '/' run_number '_' datestr(t,'yyyymmdd') '.mat'], 'out2')
                    out.STRATIGRAPHY=[];
                    out.LATERAL=[];
                    out.TIMESTAMP=[];
                    out.MISC=[];
                    out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end
            end
        end
        
        
        function [result] = usableOUT(out, interpVect)
            % Functions that transform the native CG outputs into more usable ones.
            %
            % Input argument : the out file from CG
            % Ouput argument : the tidied output
            %
            % Note 1 : This functions does not define the class that are used. If
            % needed, run the first line of the main CryoGrid program (from the
            % beggining until the time integration routine, without including it) so
            % that matlab knows them.
            %
            % Note 2 : The function is for now coded for cases without Xcess ice, with
            % a fixed grid (regarding spacing, mineral and organic content) along time.
            % See the fprintf for cases discrimintation.
            %
            % Note 3 : Arbitrary choices were made regarding the variables of interest
            % and this function can provide the base for other functions to extract and
            % present other variables.
            %
            % Note 4 : the function offer the possibility to interpolate the result at
            % a new resolution. For this, uses the input interpVect variable. It has
            % to include 2 element, the resolution of the interpolation and the depth
            % to which the interpolation is stopped, both in meters.
            %          Example : interpVect = [0.02 3];
            %
            % Author : Léo Martin, l.c.p.martin@uu.nl, October 2020, Oslo
            
            % Initiate output
            nb_dates=length(out.TIMESTAMP);
            result.TIMESTAMP=out.TIMESTAMP;
            result.T = [];
            result.waterIce =[];
            result.water = [];
            result.ice = [];
            result.air=[];
            result.snow(nb_dates).layerThick=[];
            result.snow(nb_dates).area=[];
            result.snow(nb_dates).distrThick=[];
            
            % Find number of ground modules (non snow) <------------------ Find a better way to do this
            li_snow=zeros(length(out.STRATIGRAPHY{1,1}(:,1)),1);
            for i=1:length(out.STRATIGRAPHY{1,1})
                li_snow(i)=sum(strcmp('sublimation',fieldnames(out.STRATIGRAPHY{1,1}{i,1}.STATVAR)));
            end
            nb_ground=sum(1-li_snow);
            
            % Initialize upper and lower position
            ground_names=cell(1,nb_ground);
            for i=1:nb_ground
                ground_names{i}=['ground' num2str(i)];
            end
            upperPos  =cell(nb_dates,nb_ground);
            lowerPos  =cell(nb_dates,nb_ground);
            layerThick=cell(nb_dates,nb_ground);
            area      =cell(nb_dates,nb_ground);
            
            % Browse and store dimension and snow
            for date_i=1:length(out.TIMESTAMP)
                
                STRAT_i=out.STRATIGRAPHY{1,date_i};
                
                % Store ground information
                for layer_i=length(STRAT_i)-nb_ground+1:1:length(STRAT_i)
                    upperPos{date_i,layer_i-(length(STRAT_i)-nb_ground)}  =STRAT_i{layer_i,1}.STATVAR.upperPos;
                    lowerPos{date_i,layer_i-(length(STRAT_i)-nb_ground)}  =STRAT_i{layer_i,1}.STATVAR.lowerPos;
                    layerThick{date_i,layer_i-(length(STRAT_i)-nb_ground)}=STRAT_i{layer_i,1}.STATVAR.layerThick;
                    area{date_i,layer_i-(length(STRAT_i)-nb_ground)}=STRAT_i{layer_i,1}.STATVAR.area;
                end
                
                %Store snow information assuming only one snow module  <----------------------------- Check if this assumption is always correct
                if length(STRAT_i)==nb_ground
                    result.snow(date_i).layerThick=0;
                    result.snow(date_i).area=0;
                    result.snow(date_i).distrThick=0;
                else
                    result.snow(date_i).layerThick=sum(STRAT_i{1,1}.STATVAR.layerThick);
                    result.snow(date_i).area=mean(STRAT_i{1,1}.STATVAR.area);
                    result.snow(date_i).distrThick=(result.snow(date_i).layerThick * result.snow(date_i).area)/STRAT_i{length(STRAT_i)-nb_ground+1,1}.STATVAR.area(end);
                end
            end
            
            % Check for trivial case with no subsidence and no grid modification
            if max(cell2mat(upperPos(:,1)))<= min(cell2mat(upperPos(:,1)))
                
                fprintf('usableOUT : No subsidence, simple depth processing\n')
                
                % Create z axis
                result.z.thick=vertcat(layerThick{1,:});
                result.z.TopCell=upperPos{1,1}-[0; cumsum(vertcat(layerThick{1,:}))];
                result.z.TopCell(end)=[];
                result.z.MidCell=result.z.TopCell-0.5.*vertcat(layerThick{1,:});
                A=1:nb_ground;
                B=layerThick(1,:);
                for i=1:nb_ground
                    B{i}=A(i).*B{i}./B{i};
                end
                result.z.grounds=vertcat(B{:});
                
                % Fill matrixes
                result.T = nan(length(result.z.grounds),length(result.TIMESTAMP));
                result.waterIce =result.T;
                result.water = result.T;
                result.ice = result.T;
                result.air = result.T;
                
                for date_i=1:length(out.TIMESTAMP)
                    
                    STRAT_i=out.STRATIGRAPHY{1,date_i};
                    
                    % Store ground information
                    for layer_i=length(STRAT_i)-nb_ground+1:1:length(STRAT_i)
                        % fprintf('\t%1.0f\n',layer_i)
                        result.T(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.T;
                        result.waterIce(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.waterIce;
                        result.water(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.water;
                        result.ice(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.ice;
                        result.air(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.air;
                    end
                    
                end
                
                % Fill soil data
                STRAT_i=out.STRATIGRAPHY{1,1};
                result.soil.organic=nan(length(result.z.grounds),1);
                result.soil.mineral=result.soil.organic;
                % Store ground information
                for layer_i=length(STRAT_i)-nb_ground+1:1:length(STRAT_i)
                    result.soil.organic(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),1)=STRAT_i{layer_i,1}.STATVAR.organic;
                    result.soil.mineral(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),1)=STRAT_i{layer_i,1}.STATVAR.mineral;
                end
                
                % Finalize dimensions and snow
                ground_names=cell(1,nb_ground);
                for i=1:nb_ground
                    ground_names{i}=['ground' num2str(i)];
                end
                result.dimensions.upperPos  =cell2struct(upperPos,ground_names,2);
                result.dimensions.lowerPos  =cell2struct(lowerPos,ground_names,2);
                result.dimensions.layerThick=cell2struct(layerThick,ground_names,2);
                % result.dimensions.area=cell2struct(area,ground_names,2);
                result.dimensions.area=area{1,1}(1);
                
                % Streamline dimensions
                result.dimensions.upperPos(2:end)=[];
                result.dimensions.lowerPos(2:end)=[];
                result.dimensions.layerThick(2:end)=[];
                % result.dimensions.area(2:end)=[];
                
                % Compute volumetric fractions
                [result.ice]=[result.ice]                  ./([result.z.thick]*[result.dimensions.area]);
                [result.water]=[result.water]              ./([result.z.thick]*[result.dimensions.area]);
                [result.waterIce]=[result.waterIce]        ./([result.z.thick]*[result.dimensions.area]);
                [result.air]=[result.air]                  ./([result.z.thick]*[result.dimensions.area]);
                [result.soil.organic]=[result.soil.organic]./([result.z.thick]*[result.dimensions.area]);
                [result.soil.mineral]=[result.soil.mineral]./([result.z.thick]*[result.dimensions.area]);
                
                % Interpolation
                if nargin > 1
                    % Check interpolation window
                    domainDepth=result.z.TopCell(1)-result.z.TopCell(end) + result.z.thick(end);
                    if interpVect(2)> domainDepth
                        fprintf('UsableOUT : Interpolation window deeper than modelled domain,\n            -> narrowed to model domain.\n')
                        interpVect(2)=domainDepth;
                    end
                    % Define depth vector
                    zInterp=(result.z.TopCell(1):(-1)*interpVect(1):result.z.TopCell(1)-interpVect(2))';
                    % interpolate everybody
                    goodZin=[result.z.TopCell(1);result.z.MidCell];
                    result.T=interp1(goodZin,[ result.T(1,:) ;result.T],zInterp);
                    result.water=interp1(goodZin,[ result.water(1,:) ;result.water],zInterp);
                    result.air=interp1(goodZin,[ result.air(1,:) ;result.air],zInterp);
                    result.ice=interp1(goodZin,[ result.ice(1,:) ;result.ice],zInterp);
                    result.waterIce=interp1(goodZin,[ result.waterIce(1,:) ;result.waterIce],zInterp);
                    result.soil.organic=interp1(goodZin,[ result.soil.organic(1,:) ;result.soil.organic],zInterp);
                    result.soil.mineral=interp1(goodZin,[ result.soil.mineral(1,:) ;result.soil.mineral],zInterp);
                    
                    result.z=zInterp;
                end
                
            else
                
                fprintf('usableOUT : this function does not handle simulation results with subsidence\n            Output only partially filled.\n')
                
            end
            
            
            function xls_out = write_excel(out)
                % XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
                
                xls_out = {'OUT','index',NaN,NaN;'OUT_all',1,NaN,NaN;'output_timestep',0.250000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'OUT_END',NaN,NaN,NaN};
            end
        end
    end
    
    
    
end

