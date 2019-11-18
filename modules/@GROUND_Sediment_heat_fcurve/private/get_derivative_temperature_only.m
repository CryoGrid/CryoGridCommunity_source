function divT = get_derivative_temperature_only(ground)

    %Read relevant data from ground struct
    heatFlux_ub = ground.TEMP.heatFlux_ub;    %should be generated in ground.get_boundary_condition_u
                                %or in the interaction

    Q = ground.TEMP.heatFlux_lb;       %should be generated in ground.get_boundary_condition_l
                                %or in the interaction

    T = ground.STATVAR.T;

    layerDepth = [ground.STATVAR.upperPos; ground.STATVAR.upperPos - cumsum(ground.STATVAR.layerThick)];
    midptDepth = ground.STATVAR.upperPos - ground.STATVAR.layerThick(1)/2 - ...
                cumsum(ground.STATVAR.layerThick);

    %determine bulk capacity and capacity
    c_temp = ground.STATVAR.c_eff;

    %assign thermal conductivity
    %k lives on the edges, i.e. on layerDepth
    %however, we never need the last one
    thermCond = ground.STATVAR.thermCond;


    %create output vector for temperature derivative
    divT = zeros(size(midptDepth));

    %calculate divergence
    %divT(i) = 1/c(i) (flux(i+1/2) - flux(i-1/2))/deltalayerDepth
    %where i is the position at a cell-midpoint, i.e. on midptDepth
    %and i +- 1/2 is the position at the cell-edge, i.e. on layerDepth

    %upper boundary condition
    %T(0) = TForcing
    %T(0) lives not on the midpoint of a ghost-cell,
    %but on the upper edge, hence the different deltaZ
    i = 1;
    divT(i) = (1/c_temp(i))*( ... %conductivity of cell
                            thermCond(i+1)*(T(i+1) - T(i)) / abs(midptDepth(i+1) - midptDepth(i)) - ... %flux lower edge of cell
                            heatFlux_ub) / ... %flux upper edge of cell
                    abs(layerDepth(i+1) - layerDepth(i)); %size of cell


    %derivative for soil - which difference quotient do we use here?
    for i = 2:length(midptDepth)-1
        divT(i)=(1/c_temp(i))*(...
                            thermCond(i+1) * (T(i+1) - T(i)) / abs(midptDepth(i+1) - midptDepth(i)) - ... %flux lower edge
                            thermCond(i) * (T(i) - T(i-1)) / abs(midptDepth(i) - midptDepth(i-1))) / ... %flux upper edge
                    abs(layerDepth(i+1) - layerDepth(i)); %size of cell
    end


    %lower boundary condition
    % k(N+1)(T(N+1)-T(N))/(midptDepth(N+1) - midptDepth(N)) = Q
    %thermCond(N+1) is implicitely given in Q
    i = length(midptDepth);
    divT(i) = (1/c_temp(i))*(...
                            Q - ... %flux lower edge
                            thermCond(i) * (T(i) - T(i-1)) / abs(midptDepth(i) - midptDepth(i-1))) / ... %flux upper edge
                            abs(layerDepth(i+1) - layerDepth(i)); %size of cell

    %disp(mean(divT))
end
