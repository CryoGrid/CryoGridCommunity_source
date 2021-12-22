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


classdef OUT_all_lateral_averages < matlab.mixin.Copyable
    
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
        
        
        
        %         function out = initialize_excel(out)
        %
        %         end
        
        function out = provide_PARA(out)
            % INITIALIZE_PARA  Initializes PARA structure.
            
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.time_average = [];
        end
        
        function out = provide_CONST(out)
            
        end
        
        function out = provide_STATVAR(out)
            
        end
        
 
        
        function out = finalize_init(out, tile)

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
                        res.LUT =[];  %remove look-up tables, runs out of memory otherwise
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
                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end
                    out2 = usableOUT(out);
                    save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '.mat'], 'out2')
                    out.STRATIGRAPHY=[];
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
%             % XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
%             
%             xls_out = {'OUT','index',NaN,NaN;'OUT_all',1,NaN,NaN;'output_timestep',0.250000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'OUT_END',NaN,NaN,NaN};
%         end
        
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
result.thermCond=[];
result.snow(nb_dates).layerThick=[];
result.snow(nb_dates).area=[];
result.snow(nb_dates).distrThick=[];
result.snow(nb_dates).rho=[];
result.snow(nb_dates).T=[];
result.snow(nb_dates).Qh=[];
result.snow(nb_dates).Qe=[];
result.snow(nb_dates).Lstar=[];
result.snow(nb_dates).albedo=[];
result.snow(nb_dates).thermCond=[];
result.snow(nb_dates).water=[];
result.snow(nb_dates).ice=[];
result.snow(nb_dates).air=[];
result.lateral.LAT_WATER_RESERVOIR(nb_dates).subsurface_run_off=[];
result.SEB(nb_dates).Qh=[];
result.SEB(nb_dates).Qe=[];
result.SEB(nb_dates).Qe_pot=[];
result.SEB(nb_dates).Lout=[];
result.SEB(nb_dates).Sout=[];
result.SEB(nb_dates).Lstar=[];
result.SEB(nb_dates).u_star=[];
result.runoff=nan(nb_dates,1);

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
        result.snow(date_i).layerThick=NaN;
        result.snow(date_i).area=NaN;
        result.snow(date_i).distrThick=NaN;
        result.snow(date_i).rho=NaN;
        result.snow(date_i).T=NaN;
        result.snow(date_i).Qh=NaN;
        result.snow(date_i).Qe=NaN;
        result.snow(date_i).Lstar=NaN;
        result.snow(date_i).albedo=NaN;
        result.snow(date_i).thermCond=NaN;
        result.snow(date_i).water=NaN;
        result.snow(date_i).ice=NaN;
        result.snow(date_i).air=NaN;
    else
        result.snow(date_i).layerThick=sum(STRAT_i{1,1}.STATVAR.layerThick);
        result.snow(date_i).area=mean(STRAT_i{1,1}.STATVAR.area);
        result.snow(date_i).distrThick=(result.snow(date_i).layerThick * result.snow(date_i).area)/STRAT_i{length(STRAT_i)-nb_ground+1,1}.STATVAR.area(end);
        result.snow(date_i).rho=mean(STRAT_i{1,1}.STATVAR.waterIce./STRAT_i{1,1}.STATVAR.layerThick ./STRAT_i{1,1}.STATVAR.area .*1000);
        result.snow(date_i).T=mean(STRAT_i{1,1}.STATVAR.T);
        result.snow(date_i).Qh=STRAT_i{1,1}.STATVAR.Qh;
        result.snow(date_i).Qe=STRAT_i{1,1}.STATVAR.Qe;
        result.snow(date_i).Lstar=STRAT_i{1,1}.STATVAR.Lstar;
        result.snow(date_i).albedo=STRAT_i{1,1}.STATVAR.albedo;
        result.snow(date_i).thermCond=mean(STRAT_i{1,1}.STATVAR.thermCond);
        result.snow(date_i).water=mean(STRAT_i{1,1}.STATVAR.water./STRAT_i{1,1}.STATVAR.volume);
        result.snow(date_i).ice=mean(STRAT_i{1,1}.STATVAR.ice./STRAT_i{1,1}.STATVAR.volume);
        result.snow(date_i).air=mean(STRAT_i{1,1}.STATVAR.air./STRAT_i{1,1}.STATVAR.volume);
    end
    
    % Find index of first ground
    if length(STRAT_i)==nb_ground
        index2store=1;
    else
        index2store=2;
    end
    % Store SEB info
    result.SEB(date_i).Qh=STRAT_i{index2store,1}.STATVAR.Qh;
    result.SEB(date_i).Qe=STRAT_i{index2store,1}.STATVAR.Qe;
    result.SEB(date_i).Qe_pot=STRAT_i{index2store,1}.STATVAR.Qe_pot;
    result.SEB(date_i).Lout=STRAT_i{index2store,1}.STATVAR.Lout;
    result.SEB(date_i).Sout=STRAT_i{index2store,1}.STATVAR.Sout;
    result.SEB(date_i).Lstar=STRAT_i{index2store,1}.STATVAR.Lstar;
    result.SEB(date_i).u_star=STRAT_i{index2store,1}.STATVAR.u_star;
    
    % Store scalar info from STATVAR
    if sum(double(strcmp(fieldnames(STRAT_i{index2store,1}.STATVAR), 'runoff'))) > 0
        result.runoff(date_i)=STRAT_i{index2store,1}.STATVAR.runoff;
    else
        result.runoff(date_i)=0;
    end
    
    % Store Lateral information (to be adjusted to the number of lateral module and the data they produce)
    if sum(double(strcmp(fieldnames(out.LATERAL{1,date_i}{1,1}.STATVAR), 'subsurface_run_off'))) > 0
        result.lateral.LAT_WATER_RESERVOIR(date_i).subsurface_run_off=out.LATERAL{1,date_i}{1,1}.STATVAR.subsurface_run_off;
    else
        result.lateral.LAT_WATER_RESERVOIR(date_i).subsurface_run_off=0;
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
    variable_i=1:nb_ground;
    snow_daily=layerThick(1,:);
    for i=1:nb_ground
        snow_daily{i}=variable_i(i).*snow_daily{i}./snow_daily{i};
    end
    result.z.grounds=vertcat(snow_daily{:});
    
    % Fill matrixes
    result.T = nan(length(result.z.grounds),length(result.TIMESTAMP));
    result.waterIce =result.T;
    result.water = result.T;
    result.ice = result.T;
    result.air = result.T;
    result.thermCond=result.T;
    
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
            result.thermCond(result.z.grounds==layer_i-(length(STRAT_i)-nb_ground),date_i)=STRAT_i{layer_i,1}.STATVAR.thermCond;
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
        result.thermCond=interp1(goodZin,[ result.thermCond(1,:) ;result.thermCond],zInterp);
        
        result.z=zInterp;
        
    end
    
    % Add permafrost elevation info
    [ result.frozenGround.ALT_ind, index_alongTime, result.frozenGround.ALT_val, depth_alongTime, z_alongTime ] = pfTable2(result.T,result.z.TopCell(1),result.z.MidCell);
    result.frozenGround.yearly(nb_dates).index=[];
    add=num2cell(index_alongTime);
    [result.frozenGround.yearly.index]=add{:};
    add=num2cell(depth_alongTime);
    [result.frozenGround.yearly.depth]=add{:};
    add=num2cell(z_alongTime);
    [result.frozenGround.yearly.z]=add{:};
    
    if ~isempty(out.PARA.time_average)
        % Compute daily average
        
        % Average ground data
        datetime_out=datetime(result.TIMESTAMP,'ConvertFrom','datenum');
        result.TIMESTAMP   =timeAverages(datetime_out,result.TIMESTAMP,out.PARA.time_average,'mean');
        result.T           =timeAverages(datetime_out,result.T        ,out.PARA.time_average,'mean');
        result.water       =timeAverages(datetime_out,result.water    ,out.PARA.time_average,'mean');
        result.air         =timeAverages(datetime_out,result.air      ,out.PARA.time_average,'mean');
        result.ice         =timeAverages(datetime_out,result.ice      ,out.PARA.time_average,'mean');
        result.waterIce    =timeAverages(datetime_out,result.waterIce ,out.PARA.time_average,'mean');
        result.thermCond   =timeAverages(datetime_out,result.thermCond,out.PARA.time_average,'mean');
        result.runoff      =timeAverages(datetime_out,result.runoff   ,out.PARA.time_average,'sum');
        
        % Average snow data
        snow_daily=struct2cell(result.snow);
        snow_daily=permute(snow_daily,[3 1 2]);
        where_empty=cellfun(@isempty,snow_daily); % Empty values replaced by NaN
        snow_daily(where_empty)={NaN};
        snow_daily=cell2mat(snow_daily);
        averaged_snow = timeAverages(datetime_out,snow_daily,out.PARA.time_average,'mean');
        snow_daily=num2cell(averaged_snow);
        result.snow=cell2struct(snow_daily,{'layerThick','area','distrThick','rho','T','Qh','Qe','Lstar','albedo','thermCond','water','ice','air'},2);
        
        % Average SEB data
        SEB_daily=struct2table(result.SEB);
        SEB_daily=table2array(SEB_daily)';
        averaged_SEB = timeAverages(datetime_out,SEB_daily,out.PARA.time_average,'mean');
        SEB_daily=num2cell(averaged_SEB');
        result.SEB=cell2struct(SEB_daily,{'Qh','Qe','Qe_pot','Lout','Sout','Lstar','u_star'},2);
        
        % Average Lateral reservoir data
        reservoir_daily=struct2table(result.lateral.LAT_WATER_RESERVOIR);
        reservoir_daily=table2array(reservoir_daily)';
        averaged_reservoir = timeAverages(datetime_out,reservoir_daily,out.PARA.time_average,'sum');
        reservoir_daily=num2cell(averaged_reservoir');
        result.lateral.LAT_WATER_RESERVOIR=cell2struct(reservoir_daily,{'subsurface_run_off'},2);
        
        % Average frozen ground data
        frozen_daily=struct2table(result.frozenGround.yearly);
        frozen_daily=table2array(frozen_daily)';
        averaged_frozen = timeAverages(datetime_out,frozen_daily,out.PARA.time_average,'mean');
        frozen_daily=num2cell(averaged_frozen');
        result.frozenGround.yearly=cell2struct(frozen_daily,{'index','depth','z'},2);
    end
    
else
    
    fprintf('usableOUT : this function does not handle simulation results with subsidence\n            Output only partially filled.\n')
    
end

end

function [ index, index_alongTime, ALT, depth_alongTime, z_alongTime ] = pfTable2(Tmat,z_surf,z_midCell)
% Find the permafrost table from a cryoGrid temperature matrix (Tmat, depth
% along rows and time along columns). The function can be used in "light
% mode", where just the temperature matrix is inputed, or using data about
% the model grid to produce serie of ALT and ALz

[~,nbt]=size(Tmat);

% Find index
frozen = Tmat<=0;
frozen2 = sum(frozen,2);
frozen3 = frozen2==nbt;
frozen4 = find(frozen3==1);

if isempty(frozen4) % No permafrost, trivial outputs
    index=NaN;
    ALT=NaN;
    index_alongTime=nan(nbt,1);
    depth_alongTime=nan(nbt,1);
    z_alongTime=nan(nbt,1);
    
else % There is permafrost
    
    % find index of frozen ground along time
    index=frozen4(1);
    Tmat=Tmat(1:index,:);
    unfrozen=Tmat>0; % Localize unfrozen ground
    unfrozen2=cumsum(unfrozen,'reverse'); % Sum logical indexing to have ones in the first unfrozen cell
    unfrozen3=unfrozen2==1; % Find the ones
    unfrozen4=find(unfrozen3); % Locate the index of the first unfrozen cell
    [row,col] = ind2sub(size(Tmat),unfrozen4); % convert in row col indices
    index_alongTime=nan(nbt,1); % Keep row indices
    index_alongTime(col)=row+1;
    
    if nargin==1 % Function used in light mode
        ALT=NaN;
        depth_alongTime=nan(nbt,1);
        z_alongTime=nan(nbt,1);
    else % Fuction used with elevation data
        
        % Attribute ALT
        ALT=z_surf - z_midCell(index);
        
        % find depth and z of frozen ground along time
        dealWithNaN=index_alongTime;
        dealWithNaN(isnan(dealWithNaN))=1;
        z_alongTime = z_midCell(dealWithNaN);
        depth_alongTime = z_surf - z_midCell(dealWithNaN);
        z_alongTime = z_alongTime.*(index_alongTime./index_alongTime); % Bring back the NaN where they should be
        depth_alongTime = depth_alongTime.*(index_alongTime./index_alongTime); % Bring back the NaN where they should be
        
    end
    
end

end

function [averaged_data] = timeAverages(tvec_datetime,array,period, method)
% Function that computes time averages. Time vector should be in the
% datetime format. If not the function makes the conversion.
%
% Use datetime_out=datetime(timevec,'ConvertFrom','datenum')
% to convert ISO proleptic time into datetime.
% Example of period: 'monthly', 'daily'
% Example of method: 'mean', 'sum'

% Function uses input with 1 row and many cols, so check this
[array_nl,array_nc]=size(array);
[tvec_nl,tvec_nc]=size(tvec_datetime);
if array_nl>array_nc
    array=array';
end
if tvec_nl>tvec_nc
    tvec_datetime=tvec_datetime';
end

if isdatetime(tvec_datetime)==0
    tvec_datetime=datetime(tvec_datetime,'ConvertFrom','datenum');
end

% Do the calculation
[~,nbdate]=size(tvec_datetime);
[~,nc]=size(array);
assert(nc==nbdate, 'timeAverages: dimension problem');
averaged_data = array2timetable(array','RowTimes',tvec_datetime);
averaged_data = retime(averaged_data, period, method);
averaged_data = timetable2table(averaged_data);
averaged_data(:,1)=[];
averaged_data =table2array(averaged_data)';

% Fit input format
if array_nl>array_nc
    averaged_data=averaged_data';
end

end