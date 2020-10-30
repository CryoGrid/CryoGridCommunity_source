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

% Store lateral data
if ~isempty(out.LATERAL)
    fprintf('usableOUT : out contain LATERAL data, only the subsurface_run_off value\n            of the first module of each time step will be read\n')
    result.subsurface_run_off=nan(1,length(out.LATERAL));
    for i=1:length(out.LATERAL)
        result.subsurface_run_off(i)=out.LATERAL{1,i}{1,1}.STATVAR.subsurface_run_off;
    end
end

% Infos about variables
result.varInfo(5).name=[];
result.varInfo(1).name='T';
result.varInfo(2).name='waterIce';
result.varInfo(3).name='water';
result.varInfo(4).name='ice';
result.varInfo(5).name='air';
result.varInfo(6).name='snow';
result.varInfo(7).name='z';
result.varInfo(8).name='organic';
result.varInfo(9).name='mineral';

result.varInfo(1).unit='degree C';
result.varInfo(2).unit='volumetric fraction';
result.varInfo(3).unit='volumetric fraction';
result.varInfo(4).unit='volumetric fraction';
result.varInfo(5).unit='volumetric fraction';
result.varInfo(6).unit='snow depth in meter';
result.varInfo(7).unit='meter';
result.varInfo(8).unit='volumetric fraction';
result.varInfo(9).unit='volumetric fraction';

if ~isempty(out.LATERAL)
    result.varInfo(10).name='subsurface_run_off';
    result.varInfo(10).unit='m3, cumulative';
end

end