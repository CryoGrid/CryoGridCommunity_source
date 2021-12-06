function [Water, WaterIce, snow_depth, idx] = SnowCover2(Water,WaterIce,Zs,snow_depth,snow_max,snow_fall,SnowDensity,WaterDensity)

%get current frozen snow cover dept
snow_dw = sum(Zs<=0.0); %lower most snow cell
snow_up = max(snow_dw-sum(WaterIce>0.0 & Zs<=0.0),1); %cell above last snow cell
old_snow_up = snow_up;
delta_snow_depth = snow_depth - max(-Zs(snow_up),0.0);

%snow depletion
if snow_up < snow_dw %snow cover exists
   
    %route melt water downward
    water_hold = 0.025;
    snow_rm = 0.25;
    meltwater = 0.0; %initlize melt trem
    for i=snow_up+1:snow_dw-1
        water_holding_cap = water_hold*(WaterIce(i)-Water(i));    
        %route melt water downward
        if Water(i) > water_holding_cap 
            meltwater = max(0.0, Water(i)-water_holding_cap);
            WaterIce(i) = WaterIce(i) - meltwater;
            Water(i) = Water(i) - meltwater;
            WaterIce(i+1) = WaterIce(i+1) + meltwater;
            Water(i+1) = Water(i+1) + meltwater;
                
            %remove snow cell if it contains less ice than snow_rm of inituial SWE
            if (WaterIce(i)-Water(i))<=snow_rm*SnowDensity/WaterDensity
                %treat rest as additional meltwater
                meltwater = WaterIce(i);
                WaterIce(i) = 0.0;
                Water(i) = 0.0;    
                WaterIce(i+1) = WaterIce(i+1) + meltwater;
                Water(i+1) = Water(i+1) + meltwater;
                
                %update upper snow indx and snow depth
                snow_up = snow_up+1;
                snow_depth = max(-Zs(snow_up),0.0) + delta_snow_depth;
            end
        end
    end
    
    %last snow cell
    i = snow_dw;
    %remove last snow cell and/or determine routable water
    if (WaterIce(i)-Water(i))<=snow_rm*SnowDensity/WaterDensity
        meltwater = WaterIce(i);
        WaterIce(i) = 0.0;
        Water(i) = 0.0;     
        %update upper snow indx and snow depth
        snow_up = snow_up+1;
        snow_depth = max(-Zs(snow_up),0.0) + delta_snow_depth;
    else
        meltwater = Water(i);     
        %consider this water as free routable water for upward routing
        Water(i) = Water(i)-meltwater;
        WaterIce(i) = WaterIce(i)-meltwater;
    end
    
    %route metlwater upward 
    if (meltwater>0.0) && (snow_up<snow_dw) 
        %copy water and ice arrays
        WaterIce_cp = WaterIce;
        Water_cp = Water;
        j = snow_dw;
        for k = snow_dw:-1:snow_up
            
            %step over previously removed snow cells
            while WaterIce_cp(j)==0.0 && j>=old_snow_up+1
                j = j-1;
            end
            
            %melt water routing
            if j>0 %only for safety reason
                WaterIce(k) = WaterIce_cp(j) + meltwater;
                Water(k) = Water_cp(j) + meltwater;
                meltwater = 0.0;
                
                %further route rest if cell is filled until an uptake
                %maximum
                UTmax = 0.7; %1.0 would be completely filled
                if WaterIce(k)>UTmax
                    meltwater = WaterIce(k)-UTmax;
                    WaterIce(k) = UTmax;
                    Water(k) = Water(k)-meltwater;
                end
            else %only for safety reason
                WaterIce(k) = 0.0;
                Water(k) = 0.0;
            end
            j = j-1;
        end
    end
end
%remove all water and ice rests above actual snow cover
WaterIce(1:snow_up) = 0.0;
Water(1:snow_up) = 0.0;

%snow accumulation
snow_depth = snow_depth + snow_fall*(WaterDensity/SnowDensity); %*(T0<=0.0); %increses by snow depth by snow fall
snow_depth = max(min(snow_depth,snow_max),0.0);
[~, new_snow_up] = min(abs(snow_depth- -Zs));

%build a new snow cell on top if snow depth filles a complet new cell
if new_snow_up<snow_up
    WaterIce(new_snow_up+1:snow_up) = SnowDensity/WaterDensity;
    Water(new_snow_up+1:snow_up) = 0.0;
    snow_up = new_snow_up;
end

%get new air-ground interface index
idx = snow_up+1;
end
