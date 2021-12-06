function [Water, WaterIce] = BucketScheme(Water,WaterIce,Mineral,Organic,T,Zs,dxp,dxs,PARA, rain_fall, evapotranspiration, dt)

soil_up = sum(Zs<=0.0)+1; %upper most ground cell
soil_dw = length(Zs); %lower most ground cell 

%route water downward
field_capacity_factor = 0.75; %fraction of porosity
residual_water_factor = 0.25;
%hydro_conductivity = 2e-8; %m/s

%Ice = WaterIce-Water;
Porosity = 1.0-Mineral-Organic;
%meltwater = 0.0; %initlize melt trem
%barrier = false;

% net P-ET
net_PminusET = rain_fall - evapotranspiration;

surface_runoff = 0;

Water_old = Water;

% 1. apply net P - ET to upper cells
if net_PminusET > 0 % precip input
    i = soil_up;
    excess_water = net_PminusET;
    while excess_water > 0 && T(i)>0 && i<soil_dw
        max_water = Porosity(i)*dxp(i);     % could be replaced with field capacity
        min_water = Porosity(i)*dxp(i)*residual_water_factor;         % could be replaced with residual water content
        
        actual_water = max( min_water, Water(i)*dxp(i)+excess_water );
        
        excess_water = max( 0, actual_water-max_water );
        Water(i) = min( max_water, actual_water ) / dxp(i);
        i=i+1;
        
    end
    if excess_water > 0
        surface_runoff = surface_runoff+excess_water;
    end
end
        
if net_PminusET < 0
    i = soil_up;
    deficit_water = net_PminusET;
    while deficit_water < 0 && T(i) && i<soil_dw
        max_water = Porosity(i)*dxp(i);     % could be replaced with field capacity
        min_water = Porosity(i)*dxp(i)*residual_water_factor;          % could be replaced with residual water content
        
        actual_water = min( max_water, Water(i)*dxp(i)+deficit_water );
        deficit_water = min( 0, actual_water-min_water );
        
        Water(i) = max( min_water, actual_water ) / dxp(i);
        i=i+1;
    end
    if deficit_water < 0
        warning( "water deficit could not be applied entirely");
    end
end

%2. gravity-driven infiltration, limited by Darcy-type flow
i = soil_up;


% h_upper = Water(soil_up:soil_dw-1).*dxp(soil_up:soil_dw-1);
% h_lower = -dxp(soil_up+1:soil_dw).*(1-Water(soil_up+1:soil_dw)); 
% Darcy_flow=zeros(length(Zs),1);
% Darcy_flow(soil_up:soil_dw-1)=-(h_upper-h_lower)./dxs(soil_up:soil_dw-1).*hydro_conductivity.*dt;
% Darcy_flow(Porosity==1.0 | Porosity==0.0)=0;
%assert( sum(Darcy_flow>0)==0 , "negative darcy flows occur" );

excess_water = 0;
while i < soil_dw && (T(i)>0) %barrier == false
    
%     if Porosity(i) < 1.0
%         %water_holding_cap = field_capacity_factor * Porosity(i); %this could also grid cell wise
%         %approximate max Darcy flow from cell to cell below
%         h1 = dxp(i)*Water(i); %hydrolic head in cell (Water(i) can contain water from cells above)
%         h2 = -dxp(i+1)+dxp(i+1)*Water(i+1); %negative hydrolic head in cell below
%         max_water_flow = (h1-h2)/dxs(i)*hydro_conductivity*dt; %this could also grid cell wise
%     else
%         %water_holding_cap = 1.0;
%         max_water_flow = 0.0; 
%     end 
      
%     %add rain fall in first cell if not frozen
%     if (i == soil_up) && (Ice(i) < 1e-7)
%         Water(i) = Water(i) + rain_fall/dxp(i);%vol frac
%         WaterIce(i) = WaterIce(i) + rain_fall/dxp(i);%vol frac
%     end   
    
%     %route free water downward
%     if (Water(i) > water_holding_cap) && (Ice(i) < 1e-7)   
%         meltwater = max(0.0, Water(i)-water_holding_cap); %vol frac
%         meltwater = min(meltwater,max_water_flow/dxp(i)); %vol frac
%         WaterIce(i) = WaterIce(i) - meltwater; %substract vol frac
%         Water(i) = Water(i) - meltwater; %substract vol frac
%         meltwater = meltwater*dxp(i); %in total water m³ 
%     end
    max_water = Porosity(i)*dxp(i)*field_capacity_factor;     % could be replaced with field capacity
    min_water = Porosity(i)*dxp(i)*residual_water_factor;         % could be replaced with residual water content
%    actual_water = max( min_water, Water(i).*dxp(i)-Darcy_flow(i-1)+Darcy_flow(i)+excess_water );
    actual_water = max( min_water, Water(i).*dxp(i)+excess_water );
    excess_water = max( 0, actual_water-max_water );
    Water(i) = min( max_water, actual_water ) ./ dxp(i);
    i = i+1;

%     %check for hydrological barrier 
%     if (Porosity(i+1)<1e-6) || (T(i+1)<0) 
%         barrier = true;
%     else
%         WaterIce(i+1) = WaterIce(i+1) + meltwater/dxp(i+1); %add vol frac
%         Water(i+1) = Water(i+1) + meltwater/dxp(i+1); %add vol frac
%         i = i+1;
%     end
end

i=i-1;
while excess_water > 0 && i >= soil_up
    
        max_water = Porosity(i)*dxp(i);     % could be replaced with field capacity
        min_water = Porosity(i)*dxp(i)*residual_water_factor;          % could be replaced with residual water content
        
        actual_water = max( min_water, Water(i)*dxp(i)+excess_water );
        
        excess_water = max( 0, actual_water-max_water );
        Water(i) = min( max_water, actual_water ) ./ dxp(i);
        i=i-1;
end

surface_runoff = surface_runoff + (excess_water>0)*excess_water;


% %route metlwater upward
% if meltwater>0.0
%     %copy water and ice arrays
%     WaterIce_cp = WaterIce;
%     Water_cp = Water;
%     for k = i:-1:soil_up      
%         %melt water routing
%         WaterIce(k) = WaterIce_cp(k) + meltwater/dxp(k); %add vol frac
%         Water(k) = Water_cp(k) + meltwater/dxp(k); %add vol frac
%         meltwater = 0.0;
%         
%         %further route rest if cell is filled until an uptake maximum
%         MobileWaterSpace = 1.0-Mineral(k)-Organic(k)-Ice(k);
%         UTmax = MobileWaterSpace; %pores completely filled
%         if Water(k)>UTmax
%             meltwater = Water(k)-UTmax;
%             WaterIce(k) = WaterIce(k)-meltwater;
%             Water(k) = Water(k)-meltwater;
%             meltwater = meltwater*dxp(k); %in total water m³ 
%         end        
%     end
% end

WaterIce(T>0)=Water(T>0);

% calculate storage term
WB.P=rain_fall;
WB.E = evapotranspiration;
WB.S = sum( (Water-Water_old).*dxp );
WB.R = surface_runoff;
WB.res = WB.P-WB.E-WB.R-WB.S;
if abs(WB.res)<1e-6
    fprintf("\t WB matches: P=%0.4f, ET=%0.4f, S=%0.4f, R=%0.4f\n", WB.P, WB.E,WB.S,WB.R);
else
    warning("WB mismatch!\n");
end

