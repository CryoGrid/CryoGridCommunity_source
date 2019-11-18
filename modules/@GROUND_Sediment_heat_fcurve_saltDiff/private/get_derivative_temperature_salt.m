function [divT, divsaltConc] = get_derivative_temperature_salt(ground)

    %read current state
    T = ground.STATVAR.T;
    saltConc = ground.STATVAR.saltConc;

    %read boundary conditions
    saltFlux_ub = ground.TEMP.saltFlux_ub;
    T_flux_ub = ground.TEMP.heatFlux_ub;

    saltFlux_lb = ground.TEMP.saltFlux_lb;
    T_flux_lb = ground.TEMP.heatFlux_lb;

    %read constants and grid;
    thermCond = ground.STATVAR.thermCond;
    L_f = ground.CONST.L_f;
    saltDiff = ground.STATVAR.saltDiff;

   	layerThick = ground.STATVAR.layerThick;
    midptThick = (layerThick(1:end-1)+layerThick(2:end))/2;
    %cT_delta = (ground.STATVAR.midptDepth(2:end) - ground.STATVAR.midptDepth(1:end-1));


    %heat flux divergence
    div_F_T = zeros(size(T,1),1);
    % upper boundary condition and assuming that k constant near the boundary
    div_F_T(1)= (thermCond(2).*(T(2)-T(1))./midptThick(1) ...   %Flux over first layer boundary
                - T_flux_ub)...                                 %Flux over upper boundary
                ./layerThick(1);                                %distance between boundaries

    %use FD for 2nd derivation in space
    div_F_T(2:end-1)=  (thermCond(3:end-1).*(T(3:end)-T(2:end-1))./midptThick(2:end)...       %Flux over upper layer boundary
                        - thermCond(2:end-2).*(T(2:end-1)- T(1:end-2))./midptThick(1:end-1))...   %Flux over lower layer boundary
                        ./layerThick(2:end-1);                                                  %Distance between layers

    %lower BC (dT_dt=geothermal heat flux) and assuming that k constant near the boundary
    div_F_T(end)= (T_flux_lb ....                                             %Flux over lower boundary
                    - thermCond(end-1).*(T(end)-T(end-1))./midptThick(end))...    %Flux over upper layer boundary
                    ./layerThick(end);                                          %Distance between layers


    %ion flux divergence
    div_F_saltConc=zeros(size(div_F_T));
    % upper boundary condition and assuming that k constant near the boundary
    div_F_saltConc(1)= (saltDiff(2).*(saltConc(2)-saltConc(1))./midptThick(1) ...   %Salt flux over first layer boundary
                        - saltFlux_ub)...                                           %Salt flux over upper boundary
                        ./layerThick(1);                                            %Distance between layers

    % use FD for 2nd derivation in space
    div_F_saltConc(2:end-1)=  (saltDiff(3:end-1).*(saltConc(3:end)-saltConc(2:end-1))./midptThick(2:end)...             %Salt flux over lower layer boundary
                                - saltDiff(2:end-2).*(saltConc(2:end-1)-saltConc(1:end-2))./midptThick(1:end-1))...     %Salt flux over uppper layer boundary
                                ./layerThick(2:end-1);                                                                  %Distance between layers

    % lower BC zero flux and assuming that D_C constant near the boundary
    div_F_saltConc(end)= (saltFlux_lb...                                                            %Salt flux over lower boundary
                            - saltDiff(end-1).*(saltConc(end)-saltConc(end-1))./midptThick(end))... %Salt flux over upper layer boundary
                            ./layerThick(end);                                                      %Distance between layers


    %Derivative of water content
    [dliqWater_dT, dliqWater_saltConc] = getDerivativeWaterContent(ground);
    liqWater = ground.STATVAR.liqWater;

    %put everything together

    % deriv_T_c=[ (ns./liqWater.*dliqWater_dn + 1) ./gamma .*div_F_T - Lf./gamma./liqWater .*dliqWater_dn .*div_F_ns; ...     %temperature derivative
    %         (ns./liqWater./Lf./beta .*(c./c_eff-1) .*div_F_T + 1./beta./liqWater.*div_F_ns)];                       %salt concentration derivative

    A = ground.STATVAR.c_eff;%c_eff;
    B = L_f.*dliqWater_saltConc;
    D = div_F_T;
    E = (liqWater+saltConc.*dliqWater_saltConc);
    F = saltConc.*dliqWater_dT;
    G = div_F_saltConc;

    divT = (-B.*G + D.*E)./(A.*E - B.*F);

    divsaltConc = (-F.*D + A.*G)./(A.*E - B.*F);%(-F.*D./A + G)./(E - F.*B./A);
