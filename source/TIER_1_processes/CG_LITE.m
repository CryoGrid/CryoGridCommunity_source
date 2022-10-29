%========================================================================
% CryoGrid TIER1 library class for functions related to CGLite

% S. Westermann, October 2020
%========================================================================

classdef CG_LITE < BASE
    
    methods
        function [H, dHdT] = Enthalpie(ground, T, Water, HC)
            
            rho = 1000; %[kg/m^3]
            Lsl = 334000; %[J/kg]
            L = rho*Lsl;%[J/m^3]
            
            theta = Water; %[Vol. fraction]
            
            H = (T.*HC + theta.*L); %[J/m^3]
            
            %dHdT = zeros(size(T));
            dHdT = HC.*(T~=0.0) + 1e8.*(T==0.0); %[J/m^3K]
            %dHdT = HC; %[J/m^3K]
        end
        
        function [T, fl] = EnthalpieInv(ground, H, WaterIce, HC)
            
            rho = 1000; %[kg/m^3]
            Lsl = 334000; %[J/kg]
            L = rho*Lsl;%[J/m^3]
            
            theta = WaterIce;
            %avoid zeros division
            theta(theta==0)=1e-8;
            
            fl = (H>0 & H<=L*theta).*(H./(L*theta)) + (H>L*theta);
            %c = theta.*(fl.*cl + (1-fl).*cs) + (1-theta).*cm;
            c = HC;
            
            T = (H>=(L*theta)).*((H-L*theta)./(c)) + (H<0).*(H./(c));
            
        end
        
        function HC = HeatCapacity(ground, WaterIce, Water, Mineral, Organic)
            
            cw = 4.2*10^6; %[J/m^3K] heat capacity water
            co = 2.5*10^6; %[J/m^3K]  heat capacity organic
            cm = 2.0*10^6; %[J/m^3K]  heat capacity mineral
            ca = 0.00125*10^6;%[J/m^3K]  heat capacity pore space
            ci = 1.9*10^6;%[J/m^3K]  heat capacity ice
            
            n = length(WaterIce);
            HC = zeros(n,1);
            for i=1:n
                air = 1.0 - WaterIce(i) - Mineral(i) - Organic(i);
                ice = WaterIce(i) - Water(i);
                HC(i,1) = Water(i).*cw + ice.*ci + Mineral(i).*cm + Organic(i).*co + air.*ca;
            end
            
        end
        
        function [Zn, Zs, dxn, dxs, kn, ks] = MakeGrid(ground, Zp,dxp,kp)
            
            N = length(Zp);
            Zn = Zp-dxp/2;
            Zs = Zp+dxp/2;
            
            dxn=ones(N,1);
            dxs=ones(N,1);
            for i=2:N-1
                dxn(i,1) = Zp(i) - Zp(i-1);
                dxs(i,1) = Zp(i+1) - Zp(i);
            end
            
            dxs(1,1) = Zp(2) - Zp(1);
            dxs(N,1) = dxp(N)/2;
            
            dxn(1,1) = dxp(1,1)/2;
            dxn(N,1) = Zp(N) - Zp(N-1);
            
            kn=ones(N,1);
            ks=ones(N,1);
            for i=2:N-1
                kn(i,1) = (dxp(i,1)/(2*dxn(i))*kp(i,1).^-1 + dxp(i-1,1)/(2*dxn(i))*kp(i-1).^-1).^-1;
                ks(i,1) = (dxp(i,1)/(2*dxs(i))*kp(i,1).^-1 + dxp(i+1,1)/(2*dxs(i))*kp(i+1).^-1).^-1;
            end
            kn(1,1) = kp(1);
            kn(N,1) = kp(N);
            ks(1,1) = kp(1);
            ks(N,1) = kp(N);
            
        end
        
        function [Water, WaterIce, snow_depth, idx] = SnowCover2(ground, Water,WaterIce,Zs,snow_depth,snow_max,snow_fall,SnowDensity,WaterDensity)
            
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
                meltwater = 0.0; %initlize melt term
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
        
        function x = TDMAsolver(ground, a,b,c,d)
            
            %a, b, c are the column vectors for the compressed tridiagonal matrix, d is the right vector
            n = length(b); % n is the number of rows
            x = zeros(n,1);
            % Modify the first-row coefficients
            c(1) = c(1) / b(1);    % Division by zero risk.
            d(1) = d(1) / b(1);    % Division by zero would imply a singular matrix.
            
            for i = 2:n-1
                temp = b(i) - a(i) * c(i-1);
                c(i) = c(i) / temp;
                d(i) = (d(i) - a(i) * d(i-1)) / temp;
            end
            
            d(n) = (d(n) - a(n) * d(n-1))/( b(n) - a(n) * c(n-1));
            
            % Now back substitute.
            x(n) = d(n);
            for i = n-1:-1:1
                x(i) = d(i) - c(i) * x(i + 1);
            end
            
        end
        
        function TC = ThermalConductivity(ground, WaterIce, Water, Mineral, Organic)
            
            ka = 0.025;      %air [Hillel(1982)]
            kw = 0.57;        %water [Hillel(1982)]
            ko = 0.25;        %organic [Hillel(1982)]
            km = 3.8;         %mineral [Hillel(1982)]
            %km = 2.8;         %average Granite [Bejan and Kraus (2003) in Dong et al. (2015)]
            ki = 2.2;         %ice [Hillel(1982)]
            
            n = length(WaterIce);
            TC = zeros(n,1);
            for i=1:n
                ice = WaterIce(i) - Water(i);
                air = 1.0 - WaterIce(i) - Mineral(i) - Organic(i);
                TC(i,1) = (Water(i).* kw.^0.5 + ice.* ki.^0.5 + Mineral(i).* km.^0.5 + Organic(i).* ko.^0.5 + air.* ka.^0.5).^2.0;
            end
            
        end


        
    end
end

